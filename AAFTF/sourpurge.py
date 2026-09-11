"""Run the sourmash fast matching kmer tool to look for obvious contaminants."""

import os
import shutil
import subprocess
import sys
import urllib
import uuid

from Bio import SeqIO
from packaging.version import Version

from AAFTF.resources import DB_Links
from AAFTF.utility import SafeRemove, calcN50, checkfile, execute, fastastats, get_samtools_version, printCMD, samtools_sort_cmd, samtools_view_bam_cmd, status


# logging - we may need to think about whether this has
# separate name for the different runfolder
# flake8: noqa: C901
def run(
    input,
    outfile,
    phylum,
    workdir=None,
    cpus=1,
    left=None,
    right=None,
    AAFTF_DB=None,
    sourdb=None,
    sourdb_type="gbk",
    kmer="31",
    mincovpct=5,
    taxonomy=False,
    debug=False,
    pipe=False,
    **kwargs,
):
    """Run the sourpurge routines to detect and remove contaminant contigs."""
    if not workdir:
        workdir = "aaftf-sourpurge_" + str(uuid.uuid4())[:8]
    if not os.path.exists(workdir):
        os.mkdir(workdir)

    bamthreads = 4
    if cpus < 4:
        bamthreads = 1

    # find reads
    forReads, revReads = (None,) * 2
    if left:
        forReads = os.path.abspath(left)
    if right:
        revReads = os.path.abspath(right)
    if not forReads:
        status("Unable to located FASTQ raw reads, low coverage will be skipped. Provide -l,--left (and if paired -r,--right) to enable low coverage filtering.")
        # sys.exit(1)

    # parse database locations
    if not sourdb:
        dbindex = "sourmash_genbank"
        if sourdb_type.lower() == "gtdb":
            dbindex = "sourmash_gtdb"
        elif sourdb_type.lower() == "gtdbrep" or sourdb_type.lower() == "gtdb_rep":
            dbindex = "sourmash_gtdbrep"
        elif sourdb_type.lower() == "gbk" or sourdb_type.lower() == "genbank":
            dbindex = "sourmash_gbk"
        else:
            status("Unknown sourdb_type value {:} use one of {}".format(sourdb_type, ["gtdb", "gtdbrep", "gbk"]))
            sys.exit(1)

        dburl = DB_Links[dbindex][0]["url"]
        dbfile = DB_Links[dbindex][0]["filename"]

        try:
            DB = os.environ["AAFTF_DB"]
        except KeyError:
            if AAFTF_DB:
                DB = AAFTF_DB
            else:
                status(f"$AAFTF_DB/{dbfile} not found, pass --sourdb")
                sys.exit(1)
        AAFTF_DB = DB
        SOUR = os.path.join(DB, dbfile)
        if not os.path.isfile(SOUR):
            try:
                status(f"{SOUR} sourmash database not found, downloading from {dburl} and renaming to {AAFTF_DB}/{dbfile}")
                urllib.request.urlretrieve(dburl, SOUR)
            except urllib.error.HTTPError as error:
                status(f"Error downloading from {dburl}: {error}")
                sys.exit(1)
        if not os.path.isfile(SOUR):
            status(f"{SOUR} sourmash database download of {dburl} failed. Manually download and rename to {AAFTF_DB}/{dbfile}")
            sys.exit(1)
    else:
        SOUR = os.path.abspath(sourdb)

    # hard coded tmpfile
    assembly_working = "assembly.fasta"
    blobBAM = "remapped.bam"
    shutil.copyfile(input, os.path.join(workdir, assembly_working))
    numSeqs, assemblySize = fastastats(os.path.join(workdir, assembly_working))
    status(f"Assembly is {numSeqs:,} contigs and {assemblySize:,} bp")

    # now filter for taxonomy with sourmash lca classify
    status("Running SourMash to get taxonomy classification for each contig")
    sour_sketch = os.path.basename(assembly_working) + ".sig"

    sour_compute = ["sourmash", "compute", "-k", kmer, "--scaled=1000", "--singleton", assembly_working]
    printCMD(sour_compute)
    subprocess.run(sour_compute, cwd=workdir, stderr=subprocess.DEVNULL)
    sour_classify = ["sourmash", "lca", "classify", "--db", SOUR, "--query", sour_sketch]
    printCMD(sour_classify)
    # output csv: ID,status,superkingdom,phylum,class,order,family,genus,species,strain
    Taxonomy = {}
    UniqueTax = []
    sourmashTSV = os.path.join(workdir, "sourmash.csv")
    with open(sourmashTSV, "w") as sour_out:
        for line in execute(sour_classify, workdir):
            sour_out.write(line)
            if not line or line.startswith("\n") or line.startswith("ID") or line.count(",") < 9:
                continue
            line = line.strip()
            cols = line.split(",")
            if "found" in cols[1]:
                idx = 1
                Taxonomy[cols[0]] = cols[idx + 1 :]
                taxClean = [x for x in cols[idx + 1 :] if x]
                UniqueTax.append("{:}".format(";".join(taxClean)))
            elif cols[1].strip() == "nomatch":
                idx = cols.index("nomatch")
                Taxonomy[cols[0]] = cols[idx + 1 :]
    UniqueTax = set(UniqueTax)
    status("Found {:} taxonomic classifications for contigs:\n{:}".format(len(UniqueTax), "\n".join(UniqueTax)))
    if taxonomy:
        sys.exit(1)
    Tax2Drop = []
    for k, v in Taxonomy.items():
        v = [x for x in v if x]  # remove empty items from list
        if debug:
            status(f"{k}\t{v}")
        if len(v) > 0:
            if not any(i in v for i in phylum):
                Tax2Drop.append(k)

    # drop contigs from taxonomy before calculating coverage
    status(f"Dropping {len(Tax2Drop)} contigs from taxonomy screen")
    sourTax = os.path.join(workdir, "sourmashed-tax-screen.fasta")
    with open(sourTax, "w") as sourtax_out:
        with open(os.path.join(workdir, assembly_working)) as infile:
            for record in SeqIO.parse(infile, "fasta"):
                if record.id not in Tax2Drop:
                    SeqIO.write(record, sourtax_out, "fasta")

    # only do coverage trimming if reads provided
    Contigs2Drop = []  # this will be empty if no reads given to gather by coverage
    if forReads:
        # check if BAM present, if so skip running
        if not os.path.isfile(os.path.join(workdir, blobBAM)):
            # index
            bwa_index = ["bwa", "index", os.path.basename(sourTax)]
            status("Building BWA index")
            printCMD(bwa_index)
            subprocess.run(bwa_index, cwd=workdir, stderr=subprocess.DEVNULL)
            # mapped reads to assembly using BWA
            bwa_cmd = [
                "bwa",
                "mem",
                "-t",
                str(cpus),
                os.path.basename(sourTax),  # assembly index base
                forReads,
            ]
            if revReads:
                bwa_cmd.append(revReads)

            # run BWA and pipe to samtools sort
            status("Aligning reads to assembly with BWA")
            printCMD(bwa_cmd)
            p1 = subprocess.Popen(bwa_cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            if get_samtools_version() >= Version("1.3"):
                # Modern samtools sort reads SAM directly from stdin
                sort_cmd = samtools_sort_cmd("-", os.path.join(workdir, blobBAM), bamthreads)
                printCMD(sort_cmd)
                p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stderr=subprocess.DEVNULL)
                p1.stdout.close()
                p2.communicate()
            else:
                # Older samtools: convert SAM→BAM first, then sort
                unsortBAM = os.path.join(workdir, "unsorted.bam")
                p2 = subprocess.Popen(samtools_view_bam_cmd("-", unsortBAM, bamthreads), cwd=workdir, stdin=p1.stdout, stderr=subprocess.DEVNULL)
                p1.stdout.close()
                p2.communicate()
                subprocess.run(samtools_sort_cmd(unsortBAM, os.path.join(workdir, blobBAM), bamthreads), stderr=subprocess.DEVNULL)
                SafeRemove(unsortBAM)

            subprocess.run(["samtools", "index", os.path.join(workdir, blobBAM)], stderr=subprocess.DEVNULL)

        # now calculate coverage from BAM file
        status("Calculating read coverage per contig")
        FastaBed = os.path.join(workdir, "assembly.bed")
        lengths = []
        with open(FastaBed, "w") as bedout:
            with open(sourTax) as SeqIn:
                for record in SeqIO.parse(SeqIn, "fasta"):
                    bedout.write(f"{record.id}\t{0}\t{len(record.seq)}\n")
                    lengths.append(len(record.seq))

        N50 = calcN50(lengths)
        Coverage = {}
        coverageBed = os.path.join(workdir, "coverage.bed")
        cov_cmd = ["samtools", "bedcov", os.path.basename(FastaBed), blobBAM]
        printCMD(cov_cmd)
        with open(coverageBed, "w") as bed_out:
            for line in execute(cov_cmd, workdir):
                bed_out.write(line)

                if not line or line.startswith("\n") or line.count("\t") < 3:
                    continue

                line = line.strip()
                cols = line.split("\t")
                cov = int(cols[3]) / float(cols[2])
                Coverage[cols[0]] = (int(cols[2]), cov)

        # get average coverage of N50 contigs
        n50Cov = []
        for k, v in Coverage.items():
            if debug:
                print(f"{k}; Len: {v[0]}; Cov: {v[1]:.2f}")
            if v[0] >= N50:
                n50Cov.append(v[1])
        n50AvgCov = sum(n50Cov) / len(n50Cov)
        minpct = mincovpct / 100
        # should we make this a variable? 5% was something arbitrary
        min_coverage = float(n50AvgCov * minpct)
        status(f"Average coverage for N50 contigs is {int(n50AvgCov)}X")

        # Start list of contigs to drop
        for k, v in Coverage.items():
            if v[1] <= min_coverage:
                Contigs2Drop.append(k)
        status(f"Found {len(Contigs2Drop):,} contigs with coverage less than {min_coverage:.2f}X ({mincovpct}%)")

    if debug:
        print("Contigs dropped due to coverage: {:}".format(",".join(Contigs2Drop)))
        print("Contigs dropped due to taxonomy: {:}".format(",".join(Tax2Drop)))

    DropFinal = Contigs2Drop + Tax2Drop
    DropFinal = set(DropFinal)
    status(f"Dropping {len(DropFinal):,} total contigs based on taxonomy and coverage")
    with open(outfile, "w") as outfh, open(sourTax) as seqin:
        for record in SeqIO.parse(seqin, "fasta"):
            if record.id not in DropFinal:
                SeqIO.write(record, outfh, "fasta")

    numSeqs, assemblySize = fastastats(outfile)
    status(f"Sourpurged assembly is {numSeqs:,} contigs and {assemblySize:,} bp")
    if "_" in outfile:
        nextOut = outfile.split("_")[0] + ".rmdup.fasta"
    elif "." in outfile:
        nextOut = outfile.split(".")[0] + ".rmdup.fasta"
    else:
        nextOut = outfile + ".rmdup.fasta"

    if checkfile(sourmashTSV):
        baseinput = os.path.basename(input)
        basedir = os.path.dirname(input)
        if "." in baseinput:
            baseinput = baseinput.rsplit(".", 1)[0]

        shutil.copy(sourmashTSV, os.path.join(basedir, baseinput + ".sourmash-taxonomy.csv"))

    if not debug:
        SafeRemove(workdir)

    if not pipe:
        status(f"Your next command might be:\n\tAAFTF rmdup -i {outfile} -o {nextOut}\n")
