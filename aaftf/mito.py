"""Run the Mitochondria assembly tool NOVOPlasty."""

import logging
import os
import subprocess
import sys
import uuid
from importlib.resources import as_file, files
from pathlib import Path
from typing import Any

from Bio.SeqIO.FastaIO import SimpleFastaParser

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
    memory: int = 8,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Assemble a mitochondrial genome from paired reads with NOVOPlasty.

    Two sequences steer the run:

    * the **seed** is where NOVOPlasty starts assembling; it extends the seed with matching
      reads until the genome is complete (bundled default: an *A. nidulans* cob fragment,
      ``aaftf/data/mito-seed.fasta``);
    * the **start gene** only matters afterwards: a circular assembly is rotated to begin at it
      (bundled default: a consensus of fungal cob genes, ``aaftf/data/mito-start-cob.fasta``).

    A circularized assembly is written as a single rotated ``mt`` record; otherwise the contigs
    are written as ``contig_N``.

    Args:
        left: Forward reads FASTQ.
        right: Reverse reads FASTQ.
        out: Output FASTA path.
        workdir: Working directory; a temporary one is created when None.
        minlen: Minimum expected genome size (NOVOPlasty genome range).
        maxlen: Maximum expected genome size (NOVOPlasty genome range).
        seed: FASTA of the NOVOPlasty seed; None uses the bundled *A. nidulans* cob fragment.
        starting: FASTA of the start gene a circular genome is rotated to begin at; None uses
            the bundled consensus of fungal cob genes.
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

    # NOVOPlasty seed: --seed, else the bundled default (copied into the work directory so
    # NOVOPlasty gets a real file path even from a zipped install)
    if seed:
        seed_fasta = str(Path(seed).resolve())
    else:
        seed_fasta = str(Path(workdir, "mito-seed.fasta").resolve())
        Path(seed_fasta).write_bytes((_PACKAGE_DATA / "mito-seed.fasta").read_bytes())

    # write the NOVOPlasty config from the bundled template
    placeholders = {
        "<PROJECT>": unique_id,
        "<MINLEN>": str(minlen),
        "<MAXLEN>": str(maxlen),
        "<MAXMEM>": str(memory),
        "<SEED>": seed_fasta,
        "<READLEN>": str(read_len),
        "<FORWARD>": str(Path(left).resolve()),
        "<REVERSE>": str(Path(right).resolve()),
    }
    config_text = (_PACKAGE_DATA / "novoplasty-config.txt").read_text()
    for placeholder, value in placeholders.items():
        config_text = config_text.replace(placeholder, value)
    Path(workdir, "novo-config.txt").write_text(config_text)

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
        logger.info(f"Rotating assembly to start at the start gene in {starting or 'the bundled fungal cob consensus'}")
        _orient_to_start(draft_mito, out, start_gene=starting)
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


def _orient_to_start(fasta_in: str, fasta_out: str, start_gene: str | None = None) -> None:
    """Rotate a circular MT assembly to begin at the start gene and write it as ``mt``.

    The start gene is aligned to the assembly with minimap2; with exactly one in-range hit
    the sequence is rotated (and reverse complemented for a minus-strand hit). Otherwise the
    last sequence of ``fasta_in`` is written unrotated and an error is logged.

    Args:
        fasta_in: FASTA of the circular assembly.
        fasta_out: Output FASTA path.
        start_gene: FASTA of the start gene; None uses the bundled consensus of fungal cob genes
            (``aaftf/data/mito-start-cob.fasta``).
    """
    # load sequence into dictionary
    initial_seq = ""
    # header = ''
    with open(fasta_in) as infile:
        for title, seq in SimpleFastaParser(infile):
            initial_seq = seq
            # header = title

    # without --starting, align the bundled spoa consensus of fungal cob genes; as_file gives it a
    # real path (a temporary copy only when the package is installed zipped)
    with as_file(_PACKAGE_DATA / "mito-start-cob.fasta") as default_start:
        start_file = start_gene or str(default_start)
        alignments = list(paf_hits(["minimap2", "-x", "map-ont", "-c", fasta_in, start_file]))
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
            rotated = rotated.translate(COMPLEMENT)[::-1]  # reverse complement
        with open(fasta_out, "w") as outfile:
            write_fasta(outfile, "mt", rotated)
    elif len(alignments) == 0:
        logger.error("unable to rotate because the start gene was not found in the assembly")
        with open(fasta_out, "w") as outfile:
            write_fasta(outfile, "mt", initial_seq)
    elif len(alignments) > 1:
        logger.error("unable to rotate because found multiple alignments")
        for x in alignments:
            sys.stderr.write(f"{x}\n")
        with open(fasta_out, "w") as outfile:
            write_fasta(outfile, "mt", initial_seq)
