"""Run the sourmash fast matching kmer tool to look for obvious contaminants."""

import os
import shutil
import subprocess
import sys
import urllib
import uuid
from pathlib import Path

from Bio import SeqIO

from AAFTF.resources import DB_Links
from AAFTF.utility import SafeRemove, align_to_sorted_bam, calcN50, checkfile, execute, fastastats, filter_fasta, printCMD, status


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
    if not Path(workdir).exists():
        Path(workdir).mkdir()

    bamthreads = 4
    if cpus < 4:
        bamthreads = 1

    # find reads
    forReads, revReads = (None,) * 2
    if left:
        forReads = str(Path(left).resolve())
    if right:
        revReads = str(Path(right).resolve())
    if not forReads:
        status("Unable to located FASTQ raw reads, low coverage will be skipped. Provide -l,--left (and if paired -r,--right) to enable low coverage filtering.")
        # sys.exit(1)

    # parse database locations
    if not sourdb:
        # --sourdb_type is restricted to these values by the "sourpurge" menu
        # in _menu.py (choices=["gbk", "gtdbrep", "gtdb"]).
        dbindex = {"gbk": "sourmash_gbk", "gtdbrep": "sourmash_gtdbrep", "gtdb": "sourmash_gtdb"}[sourdb_type]

        dburl = DB_Links[dbindex][0]["url"]
        dbfile = DB_Links[dbindex][0]["filename"]

        DB = os.environ.get("AAFTF_DB")
        if not DB:
            status(f"$AAFTF_DB/{dbfile} not found, pass --sourdb")
            sys.exit(1)
        SOUR = str(Path(DB, dbfile))
        if not Path(SOUR).is_file():
            try:
                status(f"{SOUR} sourmash database not found, downloading from {dburl} and renaming to {DB}/{dbfile}")
                urllib.request.urlretrieve(dburl, SOUR)
            except urllib.error.HTTPError as error:
                status(f"Error downloading from {dburl}: {error}")
                sys.exit(1)
        if not Path(SOUR).is_file():
            status(f"{SOUR} sourmash database download of {dburl} failed. Manually download and rename to {DB}/{dbfile}")
            sys.exit(1)
    else:
        SOUR = str(Path(sourdb).resolve())

    # hard coded tmpfile
    assembly_working = "assembly.fasta"
    blobBAM = "remapped.bam"
    shutil.copyfile(input, str(Path(workdir, assembly_working)))
    numSeqs, assemblySize = fastastats(str(Path(workdir, assembly_working)))
    status(f"Assembly is {numSeqs:,} contigs and {assemblySize:,} bp")

    # now filter for taxonomy with sourmash lca classify
    status("Running SourMash to get taxonomy classification for each contig")
    sour_sketch = Path(assembly_working).name + ".sig"

    sour_compute = ["sourmash", "compute", "-k", kmer, "--scaled=1000", "--singleton", assembly_working]
    printCMD(sour_compute)
    subprocess.run(sour_compute, cwd=workdir, stderr=subprocess.DEVNULL)
    sour_classify = ["sourmash", "lca", "classify", "--db", SOUR, "--query", sour_sketch]
    printCMD(sour_classify)
    # output csv: ID,status,superkingdom,phylum,class,order,family,genus,species,strain
    Taxonomy = {}
    UniqueTax = []
    sourmashTSV = str(Path(workdir, "sourmash.csv"))
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
    sourTax = str(Path(workdir, "sourmashed-tax-screen.fasta"))
    tax_drop = set(Tax2Drop)
    filter_fasta(str(Path(workdir, assembly_working)), sourTax, lambda seq_id: seq_id not in tax_drop)

    # only do coverage trimming if reads provided
    Contigs2Drop = []  # this will be empty if no reads given to gather by coverage
    if forReads:
        # check if BAM present, if so skip running
        if not Path(workdir, blobBAM).is_file():
            # index
            bwa_index = ["bwa", "index", Path(sourTax).name]
            status("Building BWA index")
            printCMD(bwa_index)
            subprocess.run(bwa_index, cwd=workdir, stderr=subprocess.DEVNULL)
            # mapped reads to assembly using BWA
            bwa_cmd = [
                "bwa",
                "mem",
                "-t",
                str(cpus),
                Path(sourTax).name,  # assembly index base
                forReads,
            ]
            if revReads:
                bwa_cmd.append(revReads)

            status("Aligning reads to assembly with BWA")
            align_to_sorted_bam(bwa_cmd, str(Path(workdir, blobBAM)), bamthreads, cwd=workdir, stderr=subprocess.DEVNULL)

            subprocess.run(["samtools", "index", str(Path(workdir, blobBAM))], stderr=subprocess.DEVNULL)

        # now calculate coverage from BAM file
        status("Calculating read coverage per contig")
        FastaBed = str(Path(workdir, "assembly.bed"))
        lengths = []
        with open(FastaBed, "w") as bedout:
            with open(sourTax) as SeqIn:
                for record in SeqIO.parse(SeqIn, "fasta"):
                    bedout.write(f"{record.id}\t{0}\t{len(record.seq)}\n")
                    lengths.append(len(record.seq))

        N50 = calcN50(lengths)
        Coverage = {}
        coverageBed = str(Path(workdir, "coverage.bed"))
        cov_cmd = ["samtools", "bedcov", Path(FastaBed).name, blobBAM]
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
    numSeqs, assemblySize = filter_fasta(sourTax, outfile, lambda seq_id: seq_id not in DropFinal)
    status(f"Sourpurged assembly is {numSeqs:,} contigs and {assemblySize:,} bp")
    if "_" in outfile:
        nextOut = outfile.split("_")[0] + ".rmdup.fasta"
    elif "." in outfile:
        nextOut = outfile.split(".")[0] + ".rmdup.fasta"
    else:
        nextOut = outfile + ".rmdup.fasta"

    if checkfile(sourmashTSV):
        baseinput = Path(input).name
        basedir = str(Path(input).parent)
        if "." in baseinput:
            baseinput = baseinput.rsplit(".", 1)[0]

        shutil.copy(sourmashTSV, str(Path(basedir, baseinput + ".sourmash-taxonomy.csv")))

    if not debug:
        SafeRemove(workdir)

    if not pipe:
        status(f"Your next command might be:\n\tAAFTF rmdup -i {outfile} -o {nextOut}\n")
