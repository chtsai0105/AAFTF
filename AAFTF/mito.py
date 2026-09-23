"""Run the Mitochondria assembly tool NOVOPlasty."""

import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser

from AAFTF.resources import Mitoseqs
from AAFTF.utility import COMPLEMENT, estimate_read_length, execute, printCMD, status, write_fasta


def run(
    left,
    right,
    out,
    workdir=None,
    minlen=10000,
    maxlen=100000,
    seed=None,
    starting=None,
    reference=None,
    memory=8,
    pipe=False,
    **kwargs,
):
    """Run the NOVOplasty tool."""
    # first check if NOVOplasty and minimap2 are installed, else exit
    programs = ["NOVOPlasty.pl", "minimap2"]
    for x in programs:
        if not shutil.which(x):
            status(f"ERROR: {x} is not installed, exiting")
            sys.exit(1)
    # first we need to generate working directory
    unique_id = str(uuid.uuid4())[:8]
    if not workdir:
        workdir = "mito_" + unique_id
    if not Path(workdir).is_dir():
        Path(workdir).mkdir(parents=True)

    # now estimate read lengths of FASTQ
    read_len = estimate_read_length(left)

    # check for seed sequence, otherwise write one
    if not seed:
        if not reference:
            seedFasta = str(Path(Path(__file__).parent, "data", "mito-seed.fasta").resolve())
        else:
            seedFasta = str(Path(reference).resolve())
    else:
        seedFasta = str(Path(seed).resolve())

    # now write the novoplasty config file
    defaultConfig = str(Path(Path(__file__).parent, "data", "novoplasty-config.txt"))
    novoConfig = str(Path(workdir, "novo-config.txt"))
    if reference:
        refgenome = str(Path(reference).resolve())
    else:
        refgenome = ""
    checkWords = ("<PROJECT>", "<MINLEN>", "<MAXLEN>", "<MAXMEM>", "<SEED>", "<READLEN>", "<FORWARD>", "<REVERSE>", "<REFERENCE>")
    repWords = (
        unique_id,  # project
        str(minlen),  # minlen
        str(maxlen),  # maxlen
        str(memory),  # maxRAM
        seedFasta,  # seed fasta seq
        str(read_len),  # read length
        str(Path(left).resolve()),  # forward read
        str(Path(right).resolve()),  # rev read
        refgenome,
    )  # ref genome file
    with open(novoConfig, "w") as outfile:
        with open(defaultConfig) as infile:
            for line in infile:
                for check, rep in zip(checkWords, repWords):
                    line = line.replace(check, rep)
                outfile.write(line)

    # now we can finally run NOVOplasty.pl
    status("De novo assembling mitochondrial genome using NOVOplasty")
    cmd = ["NOVOPlasty.pl", "-c", "novo-config.txt"]
    printCMD(cmd)
    novolog = str(Path(workdir, "novoplasty.log"))
    with open(novolog, "w") as logfile:
        p1 = subprocess.Popen(cmd, cwd=workdir, stdout=logfile, stderr=logfile)
        p1.communicate()

    # now parse the results
    draftMito = None
    circular = False
    for f in os.listdir(workdir):
        if f.startswith("Circularized_assembly_"):
            draftMito = str(Path(workdir, f))
            circular = True
            break
        if f.startswith("Contigs_1_"):
            draftMito = str(Path(workdir, f))
            break
        if f.startswith("Uncircularized_assemblies_"):
            draftMito = str(Path(workdir, f))
            break
    if draftMito is None:
        status("NOVOplasty did not produce an assembly - check log for errors")
        return
    if circular:
        status("NOVOplasty assembled complete circular genome")
        if starting:
            status(f"Rotating assembly to start with {starting}")
        else:
            status("Rotating assembly to start with Cytochrome b (cob) gene")
        _orient_to_start(draftMito, out, folder=workdir, start=starting)
    else:
        numContigs = 0
        contigLength = 0
        with open(out, "w") as outfile:
            with open(draftMito) as infile:
                for title, seq in SimpleFastaParser(infile):
                    numContigs += 1
                    contigLength += len(seq)
                    write_fasta(outfile, f"contig_{numContigs}", seq)
        # weird formatting here for PEP8
        status(f"NOVOplasty assembled {numContigs} contigs consisting of {contigLength:,} bp," + "but was unable to circularize genome")

    status(f"AAFTF mito complete: {out}")
    if not pipe:
        shutil.rmtree(workdir)


def _rev_comp(seq):
    """Reverse complement a DNA string, preserving case."""
    return seq.translate(COMPLEMENT)[::-1]


def _orient_to_start(fasta_in, fasta_out, folder=".", start=False):
    """Reorient the MT assembly based on a starting gene (if found)."""
    # if not starting, then use cytochrome oxidase (cob)
    startFile = str(Path(folder, f"{uuid.uuid4()}.fasta"))
    if not start:
        # generated as spoa consensus from select fungal cob genes
        # move this to a configurable file
        cob1 = Mitoseqs["COB1"]
        with open(startFile, "w") as outfile:
            write_fasta(outfile, "COB", cob1)
    else:
        shutil.copyfile(start, startFile)

    # load sequence into dictionary
    initial_seq = ""
    # header = ''
    with open(fasta_in) as infile:
        for title, seq in SimpleFastaParser(infile):
            initial_seq = seq
            # header = title

    alignments = []
    minimap2_cmd = ["minimap2", "-x", "map-ont", "-c", fasta_in, startFile]
    for line in execute(minimap2_cmd):
        cols = line.rstrip().split("\t")
        alignments.append(cols)
    if len(alignments) == 1:
        ref_strand = cols[4]
        ref_offset = int(cols[2])
        if ref_strand == "-":
            ref_start = int(cols[8]) + ref_offset
        else:
            ref_start = int(cols[7]) - ref_offset
        if ref_start == len(initial_seq):
            ref_start = 0
        if ref_start < 0 or ref_start > len(initial_seq):
            # A partial alignment of the seed to the contig edge can push
            # this offset out of range; Python's negative-index slicing
            # would silently wrap and produce a bogus rotation instead of
            # erroring, so treat it the same as a failed rotation.
            status(f"ERROR: unable to rotate because computed rotation offset {ref_start} is out of range for sequence of length {len(initial_seq)}\n")
            with open(fasta_out, "w") as outfile:
                write_fasta(outfile, "mt", initial_seq)
            if Path(startFile).is_file():
                Path(startFile).unlink()
            return
        rotated = initial_seq[ref_start:] + initial_seq[:ref_start]
        if ref_strand == "-":
            rotated = _rev_comp(rotated)
        with open(fasta_out, "w") as outfile:
            write_fasta(outfile, "mt", rotated)
    elif len(alignments) == 0:
        status("ERROR: unable to rotate because did " + "not find --starting sequence\n")
        with open(fasta_out, "w") as outfile:
            write_fasta(outfile, "mt", initial_seq)
    elif len(alignments) > 1:
        status("ERROR: unable to rotate because found multiple alignments\n")
        for x in alignments:
            sys.stderr.write(f"{x}\n")
        with open(fasta_out, "w") as outfile:
            write_fasta(outfile, "mt", initial_seq)
    if Path(startFile).is_file():
        Path(startFile).unlink()
