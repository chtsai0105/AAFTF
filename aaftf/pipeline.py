"""Support pipelining of AAFTF to simplify all-in-one runs.

Runs trim → filter → assemble → vecscreen → sourpurge → rmdup → polish → sort → assess.
Every step runs with its own ``AAFTF <step>`` defaults; only the options given to
``AAFTF pipeline`` (and the file names that chain the steps together) override them.
"""

import argparse as ap
import functools
import logging
from types import ModuleType
from typing import Any

import aaftf.assemble as assemble
import aaftf.assess as assess
import aaftf.filter as aaftf_filter
import aaftf.polish as polish
import aaftf.rmdup as rmdup
import aaftf.sort as aaftf_sort
import aaftf.sourpurge as sourpurge
import aaftf.trim as trim
import aaftf.vecscreen as vecscreen
from aaftf.utility import check_file

__all__ = ["run"]


logger = logging.getLogger(__name__)


def run(
    read1: str,
    basename: str,
    phylum: list[str],
    read2: str | None = None,
    cpus: int = 1,
    tmpdir: str | None = None,
    assembler_args: list[str] | None = None,
    method: str = "spades",
    memory: int | None = None,
    minlen: int = 75,
    screen_accessions: list[str] | None = None,
    screen_urls: list[str] | None = None,
    mincontiglen: int = 500,
    workdir: str | None = None,
    sourdb: str | None = None,
    mincovpct: int = 5,
    debug: bool = False,
    quiet: bool = False,
    **kwargs: Any,
) -> None:
    """Run the whole AAFTF pipeline, skipping steps whose output file already exists.

    Args:
        read1: Read 1 (forward, or single-end) raw FASTQ reads.
        basename: Prefix for every output file (``{basename}_1P.fastq.gz``, ``{basename}.final.fasta``, ...).
        phylum: Phyla whose sourmash matches are kept by ``sourpurge``.
        read2: Read 2 (reverse) raw FASTQ reads, or None for single-end data.
        cpus: Threads for every step that takes ``-c/--cpus``.
        tmpdir: Assembler temporary directory.
        assembler_args: Extra arguments passed through to the assembler.
        method: Assembler to use (``spades``, ``megahit`` or ``unicycler``).
        memory: Memory in GB for steps that take ``-m/--memory``; None keeps each step's default.
        minlen: Minimum read length after trimming.
        screen_accessions: GenBank accessions to screen out of the reads in ``filter``.
        screen_urls: URLs of sequences to screen out of the reads in ``filter``.
        mincontiglen: Minimum contig length kept by ``rmdup`` and ``sort``.
        workdir: Working directory for steps that take ``-w/--workdir``; None keeps each step's default.
        sourdb: Sourmash LCA database for ``sourpurge``.
        mincovpct: Minimum percent of N50 coverage below which ``sourpurge`` removes contigs.
        debug: Show debug output and keep temporary files in every step.
        quiet: Only show warnings and errors in every step.
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ...); ignored.

    Raises:
        RuntimeError: If a step finishes without producing its expected output file.
    """
    # passed to every step that has an option of the same name
    shared = {"cpus": cpus, "memory": memory, "workdir": workdir, "debug": debug, "quiet": quiet}
    trimmed_1, trimmed_2 = basename + "_1P.fastq.gz", (basename + "_2P.fastq.gz" if read2 else None)
    filtered_1, filtered_2 = basename + "_filtered_1.fastq.gz", (basename + "_filtered_2.fastq.gz" if read2 else None)

    _run_step(trim, "trim", trimmed_1, shared, read1=read1, read2=read2, basename=basename, minlen=minlen)
    _run_step(aaftf_filter, "filter", filtered_1, shared, read1=trimmed_1, read2=trimmed_2, basename=basename, screen_accessions=screen_accessions, screen_urls=screen_urls)

    assembly = basename + f".{method}.fasta"
    _run_step(assemble, "assemble", assembly, shared, read1=filtered_1, read2=filtered_2, out=assembly, method=method, tmpdir=tmpdir, assembler_args=assembler_args)

    vecscreen_file = basename + ".vecscreen.fasta"
    _run_step(vecscreen, "vecscreen", vecscreen_file, shared, infile=assembly, outfile=vecscreen_file)

    sourpurge_file = basename + ".sourpurge.fasta"
    _run_step(
        sourpurge,
        "sourpurge",
        sourpurge_file,
        shared,
        input=vecscreen_file,
        outfile=sourpurge_file,
        read1=filtered_1,
        read2=filtered_2,
        phylum=phylum,
        sourdb=sourdb,
        mincovpct=mincovpct,
    )

    rmdup_file = basename + ".rmdup.fasta"
    _run_step(rmdup, "rmdup", rmdup_file, shared, input=sourpurge_file, out=rmdup_file, minlen=mincontiglen)

    polish_file = basename + ".polish.fasta"
    _run_step(polish, "polish", polish_file, shared, infile=rmdup_file, outfile=polish_file, read1=filtered_1, read2=filtered_2)

    final_file = basename + ".final.fasta"
    _run_step(aaftf_sort, "sort", final_file, shared, input=polish_file, out=final_file, minlen=mincontiglen)

    assess.run(**_step_kwargs("assess", shared, input=final_file))


def _run_step(module: ModuleType, name: str, output: str, shared: dict[str, Any], **step_options: Any) -> None:
    """Run ``module.run()`` for step ``name`` unless ``output`` already exists, then check that it does.

    Args:
        module: The AAFTF subcommand module whose ``run()`` to call.
        name: Subcommand name, used to look up its CLI defaults and in log messages.
        output: File the step must produce; the step is skipped if it already exists.
        shared: Pipeline-wide options (see ``_step_kwargs``).
        **step_options: Step-specific keyword arguments that override the defaults.

    Raises:
        RuntimeError: If ``output`` is missing or empty after the step runs.
    """
    if check_file(output):
        logger.info(f"AAFTF {name} output found: {output}")
        return
    module.run(**_step_kwargs(name, shared, **step_options))
    if not check_file(output):
        raise RuntimeError(f"AAFTF {name} failed: {output} is missing or empty")


def _step_kwargs(name: str, shared: dict[str, Any], **step_options: Any) -> dict[str, Any]:
    """Return the run() keyword arguments for step ``name``: its CLI defaults, overridden by the pipeline's options.

    ``shared`` options are only passed to steps that have them, and a ``None``
    shared value (e.g. no pipeline ``--memory``) keeps the step's own default. ``pipe`` is always
    set to True so the steps don't log their "next command" hints (the pipeline runs the next step
    itself); it is a ``run()``-only parameter, not a CLI option.

    Args:
        name: Subcommand name whose CLI defaults to start from.
        shared: Pipeline-wide options (``cpus``, ``memory``, ``workdir``, ``debug``, ``quiet``).
        **step_options: Step-specific keyword arguments; these override everything else.

    Returns:
        Keyword arguments ready to pass to the step module's ``run()``.
    """
    kwargs = dict(_subcommand_defaults()[name])
    for key, value in shared.items():
        if key in kwargs and value is not None:
            kwargs[key] = value
    kwargs.update(step_options)
    kwargs["pipe"] = True
    return kwargs


@functools.cache
def _subcommand_defaults() -> dict[str, dict[str, Any]]:
    """Return ``{subcommand: {dest: default}}`` read from the AAFTF subcommand parsers."""
    from aaftf._menu import register_subcommands  # _menu imports this module, so import it when first needed

    subparsers = register_subcommands(ap.ArgumentParser())
    return {name: {action.dest: action.default for action in parser._actions if action.dest != "help"} for name, parser in subparsers.choices.items()}
