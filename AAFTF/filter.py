"""Filter runs routines to remove sequence reads matching contaminant DB.

This will match contaminant database and PhiX using read mapping kmer
tools. See resources.py for these defaults.
"""

import gzip
import os
import shutil
import subprocess
import sys
import urllib.request
import uuid
from pathlib import Path

from AAFTF.resources import Contaminant_Accessions, DB_Links, SeqDBs
from AAFTF.utility import SafeRemove, bam_read_count, countfastq, getRAM, printCMD, samtools_sort_cmd, samtools_view_bam_cmd, status


# flake8: noqa: C901
def run(
    left,
    right=None,
    workdir=None,
    cpus=1,
    AAFTF_DB=None,
    screen_accessions=None,
    screen_urls=None,
    screen_local=None,
    basename=None,
    aligner="bbduk",
    memory=None,
    debug=False,
    pipe=False,
    **kwargs,
):
    """Generic run command for this submodule for filtering reads."""
    custom_workdir = 1
    if not workdir:
        custom_workdir = 0
        workdir = "aaftf-filter_" + str(uuid.uuid4())[:8]
    if not Path(workdir).exists():
        Path(workdir).mkdir()

    # parse database locations
    DB = None
    if not AAFTF_DB:
        try:
            DB = os.environ["AAFTF_DB"]
        except KeyError:
            if AAFTF_DB:
                DB = AAFTF_DB
            else:
                pass
    else:
        DB = AAFTF_DB

    bamthreads = 4
    if cpus < 4:
        bamthreads = cpus

    earliest_file_age = -1
    contam_filenames = []
    # db of contaminant (PhiX)
    for urls in Contaminant_Accessions.values():
        for url in urls:
            acc = Path(url).name
            if DB:
                acc_file = str(Path(DB, acc))
            else:
                acc_file = str(Path(workdir, acc))
            contam_filenames.append(acc_file)
        if not Path(acc_file).exists():
            try:
                urllib.request.urlretrieve(url, acc_file)
            except Exception as e:
                status(f"error with url {url} {acc_file}: {e}")
            if earliest_file_age < 0 or earliest_file_age < Path(acc_file).stat().st_ctime:
                earliest_file_age = Path(acc_file).stat().st_ctime

    # download univec
    for url in DB_Links["UniVec"]:
        # take first file for now, could combine in future
        acc = Path(url).name
        if DB:
            acc_file = str(Path(DB, acc))
        else:
            acc_file = str(Path(workdir, acc))
        contam_filenames.append(acc_file)
        if not Path(acc_file).exists():
            urllib.request.urlretrieve(url, acc_file)
            if earliest_file_age < 0 or earliest_file_age < Path(acc_file).stat().st_ctime:
                earliest_file_age = Path(acc_file).stat().st_ctime

    if screen_accessions:
        for acc in screen_accessions:
            if DB:
                acc_file = str(Path(DB, acc + ".fna"))
                if not Path(acc_file).exists():
                    acc_file = str(Path(workdir, acc + ".fna"))
            else:
                acc_file = str(Path(workdir, acc + ".fna"))
            contam_filenames.append(acc_file)
            if not Path(acc_file).exists():
                url = SeqDBs["nucleotide"] % (acc)
                urllib.request.urlretrieve(url, acc_file)
            if earliest_file_age < 0 or earliest_file_age < Path(acc_file).stat().st_ctime:
                earliest_file_age = Path(acc_file).stat().st_ctime

    if screen_urls:
        for url in screen_urls:
            url_file = str(Path(workdir, Path(url).name))
            contam_filenames.append(url_file)
            if not Path(url_file).exists():
                urllib.request.urlretrieve(url, url_file)
            if earliest_file_age < 0 or earliest_file_age < Path(url_file).stat().st_ctime:
                earliest_file_age = Path(url_file).stat().st_ctime

    if screen_local:
        for f in screen_local:
            contam_filenames.append(str(Path(f).resolve()))

    # concat vector db

    contamdb = str(Path(workdir, "contamdb.fa"))
    filelist = "\n".join(contam_filenames)
    status(f"Generating combined contamination database {contamdb} from:\n{filelist}")
    if not Path(contamdb).exists() or (Path(contamdb).stat().st_ctime < earliest_file_age):
        with open(contamdb, "wb") as wfd:
            for fname in contam_filenames:
                if fname.endswith(".gz"):
                    with gzip.open(fname, "r") as fd:
                        shutil.copyfileobj(fd, wfd)
                else:
                    with open(fname, "rb") as fd:  # reasonably fast copy for append
                        shutil.copyfileobj(fd, wfd)

    # find reads
    forReads, revReads = (None,) * 2
    if left:
        forReads = str(Path(left).resolve())
    if right:
        revReads = str(Path(right).resolve())
    if not forReads:
        status("Must provide --left, unable to locate FASTQ reads")
        sys.exit(1)
    total = countfastq(forReads)
    if revReads:
        total = total * 2
    status(f"Loading {total:,} total reads")

    # seems like this needs to be stripping trailing extension?
    if not basename:
        if "_" in Path(forReads).name:
            basename = Path(forReads).name.split("_")[0]
        elif "." in Path(forReads).name:
            basename = Path(forReads).name.split(".")[0]
        else:
            basename = Path(forReads).name

    # logger.info('Loading {:,} FASTQ reads'.format(countfastq(forReads)))

    alignBAM = str(Path(workdir, basename + "_contam_db.bam"))
    unsorted_bam = str(Path(workdir, basename + "_contam.unsorted.bam"))
    clean_reads = basename + "_filtered"
    refmatch_bbduk = [contamdb, "phix", "artifacts", "lambda"]
    if aligner == "bbduk":
        status("Kmer filtering reads using BBDuk")
        if memory:
            MEM = f"-Xmx{memory}g"
        else:
            MEM = f"-Xmx{round(0.6 * getRAM())}g"
        leftcleanfname = f"{clean_reads}_1.fastq.gz"
        if revReads:
            # Paired mode (in=/in2=) hits a bug in this BBDuk build's
            # PairStreamer on large paired FASTQ: it silently truncates the
            # stream after a few hundred reads (with or without threading),
            # instead of erroring loudly. Feed it an INTERLEAVED single file
            # instead (bbduk's single-end reader handles the full file
            # correctly), then de-interleave the cleaned output.
            interleaved_in = str(Path(workdir, f"{basename}_ivl.fq.gz"))
            interleaved_out = str(Path(workdir, f"{basename}_ivl.clean.fq.gz"))
            shuffle_cmd = ["shuffle.sh", f"in1={forReads}", f"in2={revReads}", f"out={interleaved_in}"]
            printCMD(shuffle_cmd)
            if debug:
                subprocess.run(shuffle_cmd)
            else:
                subprocess.run(shuffle_cmd, stderr=subprocess.DEVNULL)

            cmd = ["bbduk.sh", MEM, f"t={cpus}", "hdist=1", "k=27", "overwrite=true", f"in={interleaved_in}", "interleaved=true", f"out={interleaved_out}"]
            cmd.extend(["ref={}".format(",".join(refmatch_bbduk))])
            printCMD(cmd)
            if debug:
                subprocess.run(cmd)
            else:
                subprocess.run(cmd, stderr=subprocess.DEVNULL)

            reformat_cmd = ["reformat.sh", f"in={interleaved_out}", f"out1={clean_reads}_1.fastq.gz", f"out2={clean_reads}_2.fastq.gz"]
            printCMD(reformat_cmd)
            if debug:
                subprocess.run(reformat_cmd)
            else:
                subprocess.run(reformat_cmd, stderr=subprocess.DEVNULL)
        else:
            cmd = ["bbduk.sh", MEM, f"t={cpus}", "hdist=1", "k=27", "overwrite=true"]
            cmd.extend([f"in={forReads}", f"out={clean_reads}_U.fastq.gz"])
            leftcleanfname = f"{clean_reads}_U.fastq.gz"
            cmd.extend(["ref={}".format(",".join(refmatch_bbduk))])
            # cmd.extend(['prealloc','qhdist=1'])
            printCMD(cmd)
            if debug:
                subprocess.run(cmd)
            else:
                subprocess.run(cmd, stderr=subprocess.DEVNULL)

        if not debug and not custom_workdir:
            SafeRemove(workdir)

        clean = countfastq(leftcleanfname)
        if revReads:
            clean = clean * 2  # might want to actually count - but should be always 2x
        status(f"{(total - clean):,} reads mapped to contamination database")
        status(f"{clean:,} reads unmapped and writing to file")

        if revReads:
            status("Filtering complete:\n\tFor: {:}\n\tRev: {:}".format(clean_reads + "_1.fastq.gz", clean_reads + "_2.fastq.gz"))
            if not pipe:
                status("Your next command might be:\n\tAAFTF assemble -l {:} -r {:} -c {:} -o {:}\n".format(clean_reads + "_1.fastq.gz", clean_reads + "_2.fastq.gz", cpus, basename + ".spades.fasta"))

        else:
            status("Filtering complete:\n\tSingle: {:}".format(clean_reads + "_U.fastq.gz"))
            if not pipe:
                status("Your next command might be:\n\tAAFTF assemble --merged {:} -c {:} -o {:}\n".format(clean_reads + "_U.fastq.gz", cpus, basename + ".spades.fasta"))

        return

    elif aligner == "bowtie2":
        # likely not used and less accurate than bbmap?
        if not Path(alignBAM).is_file():
            status("Aligning reads to contamination database using bowtie2")
            if not Path(contamdb + ".1.bt2").exists() or Path(contamdb + ".1.bt2").stat().st_ctime < Path(contamdb).stat().st_ctime:
                # (re)build index if no index or index is older than
                # the db
                bowtie_index = ["bowtie2-build", contamdb, contamdb]
                printCMD(bowtie_index)
                subprocess.run(bowtie_index, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

            bowtie_cmd = ["bowtie2", "-x", Path(contamdb).name, "-p", str(cpus), "--very-sensitive"]
            if forReads and revReads:
                bowtie_cmd = bowtie_cmd + ["-1", forReads, "-2", revReads]
            elif forReads:
                bowtie_cmd = bowtie_cmd + ["-U", forReads]

            # now run and write to BAM sorted
            printCMD(bowtie_cmd)

            p1 = subprocess.Popen(bowtie_cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p2 = subprocess.Popen(samtools_view_bam_cmd("-", unsorted_bam, bamthreads), cwd=workdir, stdin=p1.stdout, stderr=subprocess.DEVNULL)
            p1.stdout.close()
            p2.communicate()
            subprocess.run(samtools_sort_cmd(unsorted_bam, alignBAM, bamthreads), stderr=subprocess.DEVNULL)
            SafeRemove(unsorted_bam)

    elif aligner == "bwa":
        # likely less accurate than bbduk so may not be used
        if not Path(alignBAM).is_file():
            status("Aligning reads to contamination database using BWA")
            if not Path(contamdb + ".amb").exists() or Path(contamdb + ".amb").stat().st_ctime < Path(contamdb).stat().st_ctime:
                bwa_index = ["bwa", "index", contamdb]
                printCMD(bwa_index)
                subprocess.run(bwa_index, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

            bwa_cmd = ["bwa", "mem", "-t", str(cpus), Path(contamdb).name, forReads]
            if revReads:
                bwa_cmd.append(revReads)

            # now run and write to BAM sorted
            printCMD(bwa_cmd)
            p1 = subprocess.Popen(bwa_cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p2 = subprocess.Popen(samtools_view_bam_cmd("-", unsorted_bam, bamthreads), cwd=workdir, stdin=p1.stdout, stderr=subprocess.DEVNULL)
            p1.stdout.close()
            p2.communicate()
            subprocess.run(samtools_sort_cmd(unsorted_bam, alignBAM, bamthreads), stderr=subprocess.DEVNULL)
            SafeRemove(unsorted_bam)

    elif aligner == "minimap2":
        # likely not used but may be useful for pacbio/nanopore?
        if not Path(alignBAM).is_file():
            status("Aligning reads to contamination database using minimap2")

            minimap2_cmd = ["minimap2", "-ax", "sr", "-t", str(cpus), Path(contamdb).name, forReads]
            if revReads:
                minimap2_cmd.append(revReads)

            # now run and write to BAM sorted
            printCMD(minimap2_cmd)
            p1 = subprocess.Popen(minimap2_cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p2 = subprocess.Popen(samtools_view_bam_cmd("-", unsorted_bam, bamthreads), cwd=workdir, stdin=p1.stdout, stderr=subprocess.DEVNULL)
            p1.stdout.close()
            p2.communicate()
            subprocess.run(samtools_sort_cmd(unsorted_bam, alignBAM, bamthreads), stderr=subprocess.DEVNULL)
            SafeRemove(unsorted_bam)
    else:
        status("Must specify bowtie2, bwa, or minimap2 for filtering")

    if Path(alignBAM).is_file():
        # display mapping stats in terminal
        subprocess.run(["samtools", "index", alignBAM])
        mapped, unmapped = bam_read_count(alignBAM)
        status(f"{mapped:,} reads mapped to contamination database")
        status(f"{unmapped:,} reads unmapped and writing to file")
        # now output unmapped reads from bamfile
        # this needs to be -f 5 so unmapped-pairs
        if forReads and revReads:
            samtools_cmd = ["samtools", "fastq", "-f", "12", "-1", clean_reads + "_1.fastq.gz", "-2", clean_reads + "_2.fastq.gz", alignBAM]
        elif forReads:
            samtools_cmd = ["samtools", "fastq", "-f", "4", "-1", clean_reads + ".fastq.gz", alignBAM]
        subprocess.run(samtools_cmd, stderr=subprocess.DEVNULL)
        if not debug:
            SafeRemove(workdir)

        if revReads:
            status("Filtering complete:\n\tFor: {:}\n\tRev: {:}".format(clean_reads + "_1.fastq.gz", clean_reads + "_2.fastq.gz"))
            if not pipe:
                status("Your next command might be:\n\tAAFTF assemble -l {:} -r {:} -c {:} -o {:}\n".format(clean_reads + "_1.fastq.gz", clean_reads + "_2.fastq.gz", cpus, basename + ".spades.fasta"))
        else:
            status("Filtering complete:\n\tSingle: {:}".format(clean_reads + ".fastq.gz"))
            if not pipe:
                status("Your next command might be:\n\tAAFTF assemble -l {:} -c {:} -o {:}\n".format(clean_reads + ".fastq.gz", cpus, basename + ".spades.fasta"))
