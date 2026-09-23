"""Filter runs routines to remove sequence reads matching contaminant DB.

This will match contaminant database and PhiX using read mapping kmer
tools. See resources.py for these defaults.
"""

import logging
import sys
from pathlib import Path

from AAFTF.resources import Contaminant_Accessions, DB_Links, SeqDBs
from AAFTF.utility import aaftf_db_dir, align_to_sorted_bam, bam_read_count, basename_from_reads, cleanup_workdir, concat_files, countfastq, download_file, make_workdir, run_cmd

logger = logging.getLogger(__name__)


# flake8: noqa: C901
def run(
    left,
    right=None,
    workdir=None,
    cpus=1,
    screen_accessions=None,
    screen_urls=None,
    screen_local=None,
    basename=None,
    aligner="bbduk",
    memory=8,
    debug=False,
    pipe=False,
    **kwargs,
):
    """Generic run command for this submodule for filtering reads."""
    workdir, custom_workdir = make_workdir(workdir, "filter")
    DB = aaftf_db_dir()
    bamthreads = min(cpus, 4)

    # contaminant sequences: cached in $AAFTF_DB when set, else the workdir
    contam_filenames = []
    for url in [u for urls in Contaminant_Accessions.values() for u in urls] + DB_Links["UniVec"]:
        contam_filenames.append(download_file(url, str(Path(DB or workdir, Path(url).name))))

    for acc in screen_accessions or []:
        acc_file = str(Path(DB, acc + ".fna")) if DB else ""
        if not Path(acc_file).is_file():
            acc_file = str(Path(workdir, acc + ".fna"))
        contam_filenames.append(download_file(SeqDBs["nucleotide"] % acc, acc_file))

    for url in screen_urls or []:
        contam_filenames.append(download_file(url, str(Path(workdir, Path(url).name))))

    contam_filenames.extend(str(Path(f).resolve()) for f in screen_local or [])

    contamdb = str(Path(workdir, "contamdb.fa"))
    filelist = "\n".join(contam_filenames)
    logger.info(f"Generating combined contamination database {contamdb} from:\n{filelist}")
    newest_source = max(Path(f).stat().st_ctime for f in contam_filenames)
    if not Path(contamdb).exists() or Path(contamdb).stat().st_ctime < newest_source:
        concat_files(contam_filenames, contamdb)

    # find reads
    forReads, revReads = (None,) * 2
    if left:
        forReads = str(Path(left).resolve())
    if right:
        revReads = str(Path(right).resolve())
    if not forReads:
        logger.info("Must provide --left, unable to locate FASTQ reads")
        sys.exit(1)
    total = countfastq(forReads)
    if revReads:
        total = total * 2
    logger.info(f"Loading {total:,} total reads")

    # seems like this needs to be stripping trailing extension?
    if not basename:
        basename = basename_from_reads(forReads)

    # logger.info('Loading {:,} FASTQ reads'.format(countfastq(forReads)))

    alignBAM = str(Path(workdir, basename + "_contam_db.bam"))
    clean_reads = basename + "_filtered"
    refmatch_bbduk = [contamdb, "phix", "artifacts", "lambda"]
    if aligner == "bbduk":
        logger.info("Kmer filtering reads using BBDuk")
        MEM = f"-Xmx{memory}g"
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
            run_cmd(shuffle_cmd, debug)

            cmd = ["bbduk.sh", MEM, f"t={cpus}", "hdist=1", "k=27", "overwrite=true", f"in={interleaved_in}", "interleaved=true", f"out={interleaved_out}"]
            cmd.extend([f"ref={','.join(refmatch_bbduk)}"])
            run_cmd(cmd, debug)

            reformat_cmd = ["reformat.sh", f"in={interleaved_out}", f"out1={clean_reads}_1.fastq.gz", f"out2={clean_reads}_2.fastq.gz"]
            run_cmd(reformat_cmd, debug)
        else:
            cmd = ["bbduk.sh", MEM, f"t={cpus}", "hdist=1", "k=27", "overwrite=true"]
            cmd.extend([f"in={forReads}", f"out={clean_reads}_U.fastq.gz"])
            leftcleanfname = f"{clean_reads}_U.fastq.gz"
            cmd.extend([f"ref={','.join(refmatch_bbduk)}"])
            # cmd.extend(['prealloc','qhdist=1'])
            run_cmd(cmd, debug)

        cleanup_workdir(workdir, debug, custom_workdir)

        clean = countfastq(leftcleanfname)
        if revReads:
            clean = clean * 2  # might want to actually count - but should be always 2x
        logger.info(f"{(total - clean):,} reads mapped to contamination database")
        logger.info(f"{clean:,} reads unmapped and writing to file")

        if revReads:
            logger.info(f"Filtering complete:\nFor: {clean_reads}_1.fastq.gz\nRev: {clean_reads}_2.fastq.gz")
            if not pipe:
                logger.info(f"Your next command might be:\nAAFTF assemble -l {clean_reads}_1.fastq.gz -r {clean_reads}_2.fastq.gz -c {cpus} -o {basename}.spades.fasta")

        else:
            logger.info(f"Filtering complete:\nSingle: {clean_reads}_U.fastq.gz")
            if not pipe:
                logger.info(f"Your next command might be:\nAAFTF assemble --merged {clean_reads}_U.fastq.gz -c {cpus} -o {basename}.spades.fasta")

        return

    elif aligner == "bowtie2":
        # likely not used and less accurate than bbmap?
        if not Path(alignBAM).is_file():
            logger.info("Aligning reads to contamination database using bowtie2")
            _rebuild_index_if_stale(contamdb + ".1.bt2", contamdb, ["bowtie2-build", contamdb, contamdb], debug)

            bowtie_cmd = ["bowtie2", "-x", Path(contamdb).name, "-p", str(cpus), "--very-sensitive"]
            if forReads and revReads:
                bowtie_cmd = bowtie_cmd + ["-1", forReads, "-2", revReads]
            elif forReads:
                bowtie_cmd = bowtie_cmd + ["-U", forReads]

            align_to_sorted_bam(bowtie_cmd, alignBAM, bamthreads, cwd=workdir, debug=debug)

    elif aligner == "bwa":
        # likely less accurate than bbduk so may not be used
        if not Path(alignBAM).is_file():
            logger.info("Aligning reads to contamination database using BWA")
            _rebuild_index_if_stale(contamdb + ".amb", contamdb, ["bwa", "index", contamdb], debug)

            bwa_cmd = ["bwa", "mem", "-t", str(cpus), Path(contamdb).name, forReads]
            if revReads:
                bwa_cmd.append(revReads)

            align_to_sorted_bam(bwa_cmd, alignBAM, bamthreads, cwd=workdir, debug=debug)

    elif aligner == "minimap2":
        # likely not used but may be useful for pacbio/nanopore?
        if not Path(alignBAM).is_file():
            logger.info("Aligning reads to contamination database using minimap2")

            minimap2_cmd = ["minimap2", "-ax", "sr", "-t", str(cpus), Path(contamdb).name, forReads]
            if revReads:
                minimap2_cmd.append(revReads)

            align_to_sorted_bam(minimap2_cmd, alignBAM, bamthreads, cwd=workdir, debug=debug)
    else:
        logger.info("Must specify bowtie2, bwa, or minimap2 for filtering")

    if Path(alignBAM).is_file():
        # display mapping stats in terminal
        mapped, unmapped = bam_read_count(alignBAM)
        logger.info(f"{mapped:,} reads mapped to contamination database")
        logger.info(f"{unmapped:,} reads unmapped and writing to file")
        # now output unmapped reads from bamfile
        # this needs to be -f 5 so unmapped-pairs
        if forReads and revReads:
            samtools_cmd = ["samtools", "fastq", "-f", "12", "-1", clean_reads + "_1.fastq.gz", "-2", clean_reads + "_2.fastq.gz", alignBAM]
        elif forReads:
            samtools_cmd = ["samtools", "fastq", "-f", "4", "-1", clean_reads + ".fastq.gz", alignBAM]
        run_cmd(samtools_cmd, debug)
        cleanup_workdir(workdir, debug, custom_workdir)

        if revReads:
            logger.info(f"Filtering complete:\nFor: {clean_reads}_1.fastq.gz\nRev: {clean_reads}_2.fastq.gz")
            if not pipe:
                logger.info(f"Your next command might be:\nAAFTF assemble -l {clean_reads}_1.fastq.gz -r {clean_reads}_2.fastq.gz -c {cpus} -o {basename}.spades.fasta")
        else:
            logger.info(f"Filtering complete:\nSingle: {clean_reads}.fastq.gz")
            if not pipe:
                logger.info(f"Your next command might be:\nAAFTF assemble -l {clean_reads}.fastq.gz -c {cpus} -o {basename}.spades.fasta")


def _rebuild_index_if_stale(marker_file, contamdb, build_cmd, debug):
    """(Re)build an aligner index if its marker file is missing or older than contamdb."""
    marker = Path(marker_file)
    if not marker.exists() or marker.stat().st_ctime < Path(contamdb).stat().st_ctime:
        run_cmd(build_cmd, debug, quiet_stdout=True)
