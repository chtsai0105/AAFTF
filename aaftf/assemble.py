"""Module to run a genome assembly using defaults for Fungi.

This uses SPAdes by default but additional tools like megahit are
supported and can be added. There is some access to updating
parameters but this entire package is intended to be a general
solution for draft Illumina genome processing en masse.
"""

import logging
import os
import re
import shutil
import uuid
from pathlib import Path

from aaftf.utility import fasta_stats, run_cmd

__all__ = ["run", "run_spades", "run_dipspades", "run_megahit", "run_unicycler"]


logger = logging.getLogger(__name__)


def run(
    left,
    out,
    method="spades",
    workdir=None,
    cpus=1,
    memory="32",
    isolate=True,
    careful=True,
    assembler_args=None,
    tmpdir=None,
    right=None,
    longreads=None,
    merged=None,
    haplocontigs=False,
    debug=False,
    pipe=False,
    **kwargs,
):
    """General run command for this subcommand module where parameters are consumed."""
    if method == "spades":
        run_spades(workdir=workdir, cpus=cpus, memory=memory, isolate=isolate, careful=careful, assembler_args=assembler_args, tmpdir=tmpdir, left=left, right=right, merged=merged, out=out, debug=debug, pipe=pipe)
    elif method == "dipspades":
        run_dipspades(workdir=workdir, cpus=cpus, memory=memory, assembler_args=assembler_args, haplocontigs=haplocontigs, tmpdir=tmpdir, left=left, right=right, merged=merged, out=out, debug=debug, pipe=pipe)
    elif method == "megahit":
        run_megahit(workdir=workdir, cpus=cpus, memory=memory, assembler_args=assembler_args, tmpdir=tmpdir, left=left, right=right, out=out, debug=debug, pipe=pipe)
    elif method == "masurca":
        logger.info("Masurca assembly is not yet implemented in AAFTF")
    elif method == "nextdenovo":
        logger.info("NextDenovo assembly is not yet implemented in AAFTF")
    elif method == "unicycler":
        run_unicycler(workdir=workdir, cpus=cpus, left=left, right=right, longreads=longreads, merged=merged, out=out, debug=debug, pipe=pipe)
    else:
        logger.info(f"Unknown assembler method {method}")


def run_spades(workdir=None, cpus=1, memory="32", isolate=True, careful=True, assembler_args=None, tmpdir=None, left=None, right=None, merged=None, out=None, debug=False, pipe=False, **kwargs):
    """Run SPAdes assembhler."""
    if not workdir:
        workdir = "spades_" + str(uuid.uuid4())[:8]

    runcmd = ["spades.py", "--threads", str(cpus), "--mem", memory, "-o", workdir]

    if isolate:
        runcmd.extend(["--isolate"])
    elif careful:
        runcmd.extend(["--careful"])

    if assembler_args:
        runcmd.extend(assembler_args)

    if "--meta" not in runcmd:
        runcmd.extend(["--cov-cutoff", "auto"])

    if tmpdir:
        runcmd.extend(["--tmp-dir", tmpdir])

    forward_reads, reverse_reads = _resolve_reads(left, right)

    if not reverse_reads:
        runcmd.extend(["--s1", forward_reads])
        if merged:
            runcmd.extend(["--s2", merged])
    else:
        runcmd.extend(["--pe1-1", forward_reads, "--pe1-2", reverse_reads])
        if merged:
            runcmd.extend(["--s1", merged])

    # this basically overrides everything above and only runs --restart-from option
    if Path(workdir).is_dir():
        runcmd = ["spades.py", "-o", workdir, "--threads", str(cpus), "--mem", memory, "--restart-from last"]

    logger.info("Assembling FASTQ data using Spades")
    run_cmd(runcmd, debug, quiet_stdout=True)

    final_out = _derive_final_out(out, forward_reads, ".spades.fasta")
    _finish_assembly(Path(workdir, "scaffolds.fasta"), final_out, "Spades", cpus, pipe)


def run_dipspades(workdir=None, cpus=1, memory="32", assembler_args=None, haplocontigs=False, tmpdir=None, left=None, right=None, merged=None, out=None, debug=False, pipe=False, **kwargs):
    """Run dipSPAdes for diploid assembly support, only on older version of SPAdes."""
    if not workdir:
        workdir = "dipspades_" + str(os.getpid())

    runcmd = ["dipspades.py", "--threads", str(cpus), "--cov-cutoff", "auto", "--mem", memory, "-o", workdir]

    if assembler_args:
        runcmd.extend(assembler_args)

    if haplocontigs:
        runcmd.extend(["--hap", haplocontigs])

    if tmpdir:
        runcmd.extend(["--tmp-dir", tmpdir])

    forward_reads, reverse_reads = _resolve_reads(left, right)

    if not reverse_reads:
        runcmd.extend(["-s", forward_reads])
    else:
        runcmd.extend(["--pe1-1", forward_reads, "--pe1-2", reverse_reads])
        if merged:
            runcmd.extend(["-s", merged])

    # this basically overrides everything above and only runs --restart-from option
    if Path(workdir).is_dir():
        runcmd = ["dipspades.py", "-o", workdir, "--continue"]

    logger.info("Assembling FASTQ data using Spades")
    run_cmd(runcmd, debug, quiet_stdout=True)

    final_out = _derive_final_out(out, forward_reads, ".dipspades.fasta")
    prefix = Path(final_out).name.removesuffix(".dipspades.fasta")

    if Path(workdir, "consensus_contigs.fasta").is_file():
        shutil.copyfile(str(Path(workdir, "dipspades", "paired_consensus_contigs.fasta")), prefix + ".dipspades_consensus_paired.fasta")
        shutil.copyfile(str(Path(workdir, "dipspades", "paired_consensus_contigs.fasta")), prefix + ".dipspades_consensus_unpaired.fasta")
        logger.info("Dipspades assembly copied over: {:}".format(prefix + ".dipspades_consensus_unpaired.fasta"), prefix + ".dipspades_consensus_paired.fasta")

    missing_msg = "Spades assembly output missing -- check Dipspades logfile in {:}.".format(str(Path(workdir, "dipspades", "dipspades.log")))
    _finish_assembly(Path(workdir, "consensus_contigs.fasta"), final_out, "Dipspades", cpus, pipe, missing_msg=missing_msg)


def run_megahit(workdir=None, cpus=1, memory=None, assembler_args=None, tmpdir=None, left=None, right=None, out=None, debug=False, pipe=False, **kwargs):
    """Run megahit assembler. This is faster but maybe less accurate."""
    if not workdir:
        workdir = "megahit_" + str(os.getpid())

    runcmd = ["megahit", "-t", str(cpus), "-o", workdir]

    if assembler_args:
        runcmd.extend(assembler_args)

    if memory:
        runcmd.extend(["--memory", memory])

    if tmpdir:
        runcmd.extend(["--tmp-dir", tmpdir])

    forward_reads, reverse_reads = _resolve_reads(left, right)

    if not reverse_reads:
        runcmd.extend(["-r", forward_reads])
    else:
        runcmd.extend(["-1", forward_reads, "-2", reverse_reads])

    if Path(workdir).is_dir():
        logger.info(f"Cannot re-run with existing folder {workdir}")

    logger.info("Assembling FASTQ data using megahit")
    run_cmd(runcmd, debug, quiet_stdout=True)

    final_out = _derive_final_out(out, forward_reads, ".megahit.fasta")
    _finish_assembly(Path(workdir, "final.contigs.fa"), final_out, "Megahit", cpus, pipe)


def run_unicycler(workdir=None, cpus=1, left=None, right=None, longreads=None, merged=None, out=None, debug=False, pipe=False, **kwargs):
    """Run Unicycler assembhler."""
    if not workdir:
        workdir = "unicycler_" + str(uuid.uuid4())[:8]

    runcmd = ["unicycler", "--threads", str(cpus), "-o", workdir]

    # if memory:
    #    runcmd.extend(['--spades_options', f'-m {memory}'])

    forward_reads, reverse_reads = _resolve_reads(left, right)

    if longreads:
        runcmd.extend(["--long", longreads])

    if not reverse_reads:
        runcmd.extend(["--unpaired", forward_reads])
    elif merged:
        runcmd.extend(["--unpaired", merged])
    else:
        runcmd.extend(["--short1", forward_reads, "--short2", reverse_reads])
        if merged:
            runcmd.extend(["--unpaired", merged])

    # not supporting restarting a run
    # this basically overrides everything above and only runs --restart-from option
    #    if Path(workdir).is_dir():
    #    runcmd = ['unicycler', '-o', workdir,
    #            '--threads', str(cpus),
    #            '--mem', memory,
    #            '--restart-from last']

    logger.info("Assembling FASTQ data using Unicycler")
    run_cmd(runcmd, debug, quiet_stdout=True)

    final_out = _derive_final_out(out, forward_reads, ".unicycler.fasta")
    _finish_assembly(Path(workdir, "assembly.fasta"), final_out, "Unicycler", cpus, pipe)


def _resolve_reads(left, right):
    """Resolve absolute paths for forward/reverse reads; exit if forward reads are missing."""
    forward_reads = str(Path(left).resolve()) if left else None
    reverse_reads = str(Path(right).resolve()) if right else None
    if not forward_reads:
        raise ValueError("Unable to locate FASTQ raw reads, provide --left")
    return forward_reads, reverse_reads


def _derive_final_out(out, forward_reads, suffix):
    """Derive the assembly output FASTA filename from --out, or from the input read filename."""
    if out:
        return out
    prefix = Path(forward_reads).name
    m = re.search(r"(\S+)\.(fastq|fq)(\.\S+)?", prefix)
    if m:
        prefix = m.group(1)
    return prefix + suffix


def _finish_assembly(src, final_out, tool_name, cpus, pipe, missing_msg=None):
    """Copy the assembler's raw output to finalOut, report stats, and print the next-step hint."""
    if Path(src).is_file():
        shutil.copyfile(str(src), final_out)
        logger.info(f"{tool_name} assembly finished: {final_out}")
        num_seqs, assembly_size = fasta_stats(final_out)
        logger.info(f"Assembly is {num_seqs:,} scaffolds and {assembly_size:,} bp")
    else:
        logger.info(missing_msg or f"{tool_name} assembly output missing -- check {tool_name} logfile.")

    if not pipe:
        logger.info(f"Your next command might be:\nAAFTF vecscreen -i {final_out} -c {cpus}")
