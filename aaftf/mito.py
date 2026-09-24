"""Run the Mitochondria assembly tool NOVOPlasty."""

import logging
import os
import shutil
import subprocess
import sys
import uuid
from importlib.resources import files
from pathlib import Path
from typing import Any

from Bio.SeqIO.FastaIO import SimpleFastaParser

from aaftf.resources import MITO_SEQS
from aaftf.utility import COMPLEMENT, cleanup_workdir, estimate_read_length, make_workdir, paf_hits, print_cmd, require_tools, write_fasta

__all__ = ["run"]


logger = logging.getLogger(__name__)

_PACKAGE_DATA = files("aaftf") / "data"


def run(
    left: str,
    right: str,
    out: str = "mito.fasta",
    workdir: str | None = None,
    minlen: int = 10000,
    maxlen: int = 100000,
    seed: str | None = None,
    starting: str | None = None,
    reference: str | None = None,
    memory: int = 8,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Assemble a mitochondrial genome from paired reads with NOVOPlasty.

    Writes a NOVOPlasty config from the bundled template and runs it in the work directory. A
    circularized assembly is rotated to start at ``starting`` (default: the cob gene); otherwise
    the contigs are written as ``contig_N``. Returns early, writing nothing, if NOVOPlasty
    produced no assembly.

    Args:
        left: Forward reads FASTQ.
        right: Reverse reads FASTQ.
        out: Output FASTA path.
        workdir: Working directory; a temporary one is created when None.
        minlen: Minimum expected genome size (NOVOPlasty genome range).
        maxlen: Maximum expected genome size (NOVOPlasty genome range).
        seed: Seed FASTA; falls back to ``reference``, then the bundled seed.
        starting: FASTA of the sequence to rotate a circular genome to start at.
        reference: Reference genome FASTA passed to NOVOPlasty.
        memory: Max memory in GB for NOVOPlasty.
        debug: Keep the work directory.
        pipe: Unused; accepted for pipeline consistency.
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``, ...); ignored.

    Raises:
        RuntimeError: If NOVOPlasty produces no assembly.
    """
    require_tools(["NOVOPlasty.pl", "minimap2"])
    # first we need to generate working directory
    unique_id = str(uuid.uuid4())[:8]
    workdir, custom_workdir = make_workdir(workdir, "mito")

    # now estimate read lengths of FASTQ
    read_len = estimate_read_length(left)

    # seed sequence: --seed, else --reference, else the bundled default (copied out of the
    # package so NOVOPlasty gets a real file path even from a zipped install)
    if seed:
        seed_fasta = str(Path(seed).resolve())
    elif reference:
        seed_fasta = str(Path(reference).resolve())
    else:
        seed_fasta = str(Path(workdir, "mito-seed.fasta").resolve())
        Path(seed_fasta).write_bytes((_PACKAGE_DATA / "mito-seed.fasta").read_bytes())

    # now write the novoplasty config file from the bundled template
    novo_config = str(Path(workdir, "novo-config.txt"))
    if reference:
        refgenome = str(Path(reference).resolve())
    else:
        refgenome = ""
    check_words = ("<PROJECT>", "<MINLEN>", "<MAXLEN>", "<MAXMEM>", "<SEED>", "<READLEN>", "<FORWARD>", "<REVERSE>", "<REFERENCE>")
    rep_words = (
        unique_id,  # project
        str(minlen),  # minlen
        str(maxlen),  # maxlen
        str(memory),  # maxRAM
        seed_fasta,  # seed fasta seq
        str(read_len),  # read length
        str(Path(left).resolve()),  # forward read
        str(Path(right).resolve()),  # rev read
        refgenome,
    )  # ref genome file
    config_text = (_PACKAGE_DATA / "novoplasty-config.txt").read_text()
    for check, rep in zip(check_words, rep_words):
        config_text = config_text.replace(check, rep)
    Path(novo_config).write_text(config_text)

    # now we can finally run NOVOplasty.pl
    logger.info("De novo assembling mitochondrial genome using NOVOplasty")
    cmd = ["NOVOPlasty.pl", "-c", "novo-config.txt"]
    print_cmd(cmd)
    novolog = str(Path(workdir, "novoplasty.log"))
    with open(novolog, "w") as logfile:
        p1 = subprocess.Popen(cmd, cwd=workdir, stdout=logfile, stderr=logfile)
        p1.communicate()

    # now parse the results, preferring a circular assembly over partial ones
    outputs = sorted(os.listdir(workdir))
    draft_mito = None
    for output_prefix in ("Circularized_assembly_", "Contigs_1_", "Uncircularized_assemblies_"):
        draft_mito = next((str(Path(workdir, f)) for f in outputs if f.startswith(output_prefix)), None)
        if draft_mito:
            break
    if draft_mito is None:
        raise RuntimeError(f"NOVOplasty did not produce an assembly - check {novolog}")
    circular = Path(draft_mito).name.startswith("Circularized_assembly_")
    if circular:
        logger.info("NOVOplasty assembled complete circular genome")
        if starting:
            logger.info(f"Rotating assembly to start with {starting}")
        else:
            logger.info("Rotating assembly to start with Cytochrome b (cob) gene")
        _orient_to_start(draft_mito, out, folder=workdir, start=starting)
    else:
        num_contigs = 0
        contig_length = 0
        with open(out, "w") as outfile:
            with open(draft_mito) as infile:
                for title, seq in SimpleFastaParser(infile):
                    num_contigs += 1
                    contig_length += len(seq)
                    write_fasta(outfile, f"contig_{num_contigs}", seq)
        # weird formatting here for PEP8
        logger.info(f"NOVOplasty assembled {num_contigs} contigs consisting of {contig_length:,} bp," + "but was unable to circularize genome")

    logger.info(f"AAFTF mito complete: {out}")
    cleanup_workdir(workdir, debug, custom_workdir)


def _orient_to_start(fasta_in: str, fasta_out: str, folder: str = ".", start: str | None = None) -> None:
    """Rotate a circular MT assembly to begin at a starting gene and write it as ``mt``.

    The start sequence is aligned to the assembly with minimap2; with exactly one in-range hit
    the sequence is rotated (and reverse complemented for a minus-strand hit). Otherwise the
    last sequence of ``fasta_in`` is written unrotated and an error is logged.

    Args:
        fasta_in: FASTA of the circular assembly.
        fasta_out: Output FASTA path.
        folder: Directory for the temporary start-sequence FASTA.
        start: FASTA file of the start sequence; None uses the built-in COB consensus.
    """
    # if not starting, then use cytochrome oxidase (cob)
    start_file = str(Path(folder, f"{uuid.uuid4()}.fasta"))
    if not start:
        # generated as spoa consensus from select fungal cob genes
        # move this to a configurable file
        cob1 = MITO_SEQS["COB1"]
        with open(start_file, "w") as outfile:
            write_fasta(outfile, "COB", cob1)
    else:
        shutil.copyfile(start, start_file)

    # load sequence into dictionary
    initial_seq = ""
    # header = ''
    with open(fasta_in) as infile:
        for title, seq in SimpleFastaParser(infile):
            initial_seq = seq
            # header = title

    minimap2_cmd = ["minimap2", "-x", "map-ont", "-c", fasta_in, start_file]
    try:
        alignments = list(paf_hits(minimap2_cmd))
    finally:
        Path(start_file).unlink(missing_ok=True)
    if len(alignments) == 1:
        hit = alignments[0]
        ref_strand = hit.strand
        ref_offset = hit.query_start
        if ref_strand == "-":
            ref_start = hit.target_end + ref_offset
        else:
            ref_start = hit.target_start - ref_offset
        if ref_start == len(initial_seq):
            ref_start = 0
        if ref_start < 0 or ref_start > len(initial_seq):
            # A partial alignment of the seed to the contig edge can push
            # this offset out of range; Python's negative-index slicing
            # would silently wrap and produce a bogus rotation instead of
            # erroring, so treat it the same as a failed rotation.
            logger.error(f"unable to rotate because computed rotation offset {ref_start} is out of range for sequence of length {len(initial_seq)}")
            with open(fasta_out, "w") as outfile:
                write_fasta(outfile, "mt", initial_seq)
            return
        rotated = initial_seq[ref_start:] + initial_seq[:ref_start]
        if ref_strand == "-":
            rotated = _rev_comp(rotated)
        with open(fasta_out, "w") as outfile:
            write_fasta(outfile, "mt", rotated)
    elif len(alignments) == 0:
        logger.error("unable to rotate because did " + "not find --starting sequence")
        with open(fasta_out, "w") as outfile:
            write_fasta(outfile, "mt", initial_seq)
    elif len(alignments) > 1:
        logger.error("unable to rotate because found multiple alignments")
        for x in alignments:
            sys.stderr.write(f"{x}\n")
        with open(fasta_out, "w") as outfile:
            write_fasta(outfile, "mt", initial_seq)


def _rev_comp(seq: str) -> str:
    """Reverse complement a DNA string, preserving case.

    Args:
        seq: DNA sequence.

    Returns:
        The reverse complement.
    """
    return seq.translate(COMPLEMENT)[::-1]
