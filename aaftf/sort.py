"""This module sorts FASTA sequences by size and renames headers."""

import logging

from Bio.SeqIO.FastaIO import SimpleFastaParser

from aaftf.utility import write_fasta

__all__ = ["run"]


logger = logging.getLogger(__name__)


def run(input, out, minlen=0, name="scaffold", **kwargs):
    """Sort contig/scaffold file longest to shortest and rename."""
    logger.info("Sorting sequences by length longest --> shortest")
    all_seqs = {}
    with open(input) as fasta_in:
        for header, seq in SimpleFastaParser(fasta_in):
            if header not in all_seqs:
                if len(seq) >= minlen:
                    all_seqs[header] = seq
    sorted_seqs = sorted(all_seqs.items(), key=lambda item: len(item[1]), reverse=True)
    with open(out, "w") as fasta_out:
        for i, (header, seq) in enumerate(sorted_seqs):
            write_fasta(fasta_out, f"{name}_{i + 1}", seq)

    logger.info(f"Output written to: {out}")
