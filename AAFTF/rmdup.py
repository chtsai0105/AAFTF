"""RMDup module step removes smaller contigs redudnant with larger ones.

This uses minimap to map small contigs against the database of contigs in an
assembly and removes those which are redundant.
"""

import operator
import os
import sys
import uuid
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser

from AAFTF.utility import SafeRemove, calc_nx, execute, status, write_fasta


def run(
    input,
    out,
    workdir=None,
    cpus=1,
    percent_id=95,
    percent_cov=95,
    minlen=500,
    exhaustive=False,
    debug=False,
    pipe=False,
    **kwargs,
):
    """Run routines to identify and remove duplicate contigs."""

    def generateFastas(fasta, pref, query, reference):
        qfile = str(Path(workdir, pref + "query.fasta"))
        rfile = str(Path(workdir, pref + "reference.fasta"))
        with open(qfile, "w") as qout:
            with open(rfile, "w") as rout:
                with open(fasta) as infile:
                    for Header, Seq in SimpleFastaParser(infile):
                        if Header in query:
                            write_fasta(qout, Header, Seq)
                        elif Header in reference:
                            write_fasta(rout, Header, Seq)
        return qfile, rfile

    def runMinimap2(query, reference, name):
        """Run minimap2 for matching contigs."""
        garbage = False  # assume this is a good contig
        for line in execute(["minimap2", "-t", str(cpus), "-x", "asm5", "-N5", reference, query], "."):
            qID, qLen, qStart, qEnd, strand, tID, tLen, tStart, tEnd, matches, alnLen, mapQ = line.split("\t")[:12]
            pident = float(matches) / int(alnLen) * 100
            cov = float(alnLen) / int(qLen) * 100
            if debug:
                print(f"\tquery={qID} hit={tID} pident={pident:.2f} coverage={cov:.2f}")

            if pident > percent_id and cov > percent_cov:
                if debug:
                    print(f"{name} duplicated: {pident:.0f}% identity over {cov:.0f}% of the contig. length={qLen}")
                garbage = True
                break
        return garbage  # false is good, true is repeat

    # start here -- functions nested so they can inherit the arguments
    custom_workdir = 1
    if not workdir:
        custom_workdir = 0
        workdir = "aaftf-rmdup_" + str(uuid.uuid4())[:8]
    if not Path(workdir).exists():
        Path(workdir).mkdir()

    if debug:
        status(f"input={input} out={out} workdir={workdir} cpus={cpus} percent_id={percent_id} " f"percent_cov={percent_cov} minlen={minlen} exhaustive={exhaustive} pipe={pipe}")
    status("Looping through assembly shortest --> longest searching for duplicated contigs using minimap2")
    fasta_lengths = []
    AllSeqs = {}
    with open(input) as infile:
        for Header, Seq in SimpleFastaParser(infile):
            fasta_lengths.append(len(Seq))
            AllSeqs.setdefault(Header, len(Seq))
    n50, _ = calc_nx(fasta_lengths, 0.5)
    n75, _ = calc_nx(fasta_lengths, 0.75)
    status(f"Assembly is {len(fasta_lengths):,} contigs; {sum(fasta_lengths):,} bp; N50 is {n50:,} bp; N75 is {n75:,} bp")

    # get list of tuples of sequences sorted by size (shortest --> longest)
    sortSeqs = sorted(AllSeqs.items(), key=operator.itemgetter(1), reverse=False)
    if exhaustive:
        n75 = sortSeqs[-1][1]
    those2check = [x for x in sortSeqs if x[1] < n75]
    status(f"Will check {len(those2check):,} contigs for duplication --> those that are < {n75:,} && > {minlen:,}")
    # loop through sorted list of tuples
    ignore = []
    for i, x in enumerate(sortSeqs):
        sys.stdout.flush()
        if x[1] < minlen:
            ignore.append(x[0])
            continue
        if x[1] > n75:
            sys.stdout.flush()
            sys.stdout.write("\n")
            break
        if debug:
            status(f"Working on {x[0]} len={x[1]} remove_tally={len(ignore)}")
        else:
            text = f"\rProgress: {i} of {len(those2check)}; remove tally={len(ignore):,}; current={x[0]}; length={x[1]}     "
            sys.stdout.write(text)
        # generate input files for minimap2
        theRest = [i[0] for i in sortSeqs[i + 1 :]]
        pid = str(os.getpid())
        qfile, rfile = generateFastas(input, pid, x[0], theRest)
        # run minimap2
        result = runMinimap2(qfile, rfile, x[0])
        if result:
            ignore.append(x[0])

    ignore = set(ignore)
    numSeqs = assemblySize = 0
    with open(out, "w") as clean_out:
        with open(input) as infile:
            for Header, Seq in SimpleFastaParser(infile):
                if Header not in ignore:
                    write_fasta(clean_out, Header, Seq)
                    numSeqs += 1
                    assemblySize += len(Seq)
    status(f"Cleaned assembly is {numSeqs:,} contigs and {assemblySize:,} bp")
    if "_" in out:
        nextOut = out.split("_")[0] + ".polish.fasta"
    elif "." in out:
        nextOut = out.split(".")[0] + ".polish.fasta"
    else:
        nextOut = out + ".polish.fasta"

    if not pipe:
        status(f"Your next command might be:\n\tAAFTF polish -i {out} -l PE_R1.fastq.gz -r PE_R2.fastq.gz -o {nextOut}\n")

    if not debug and not custom_workdir:
        SafeRemove(workdir)
