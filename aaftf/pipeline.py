"""Support pipelining of AAFTF to simplify all-in-one runs.

Runs trim → filter → assemble → vecscreen → sourpurge → rmdup → polish → sort → assess.
Every step runs with its own ``AAFTF <step>`` defaults; only the options given to
``AAFTF pipeline`` (and the file names that chain the steps together) override them.
"""

import argparse as ap
import functools
import logging

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
    left,
    basename,
    phylum,
    right=None,
    cpus=1,
    tmpdir=None,
    assembler_args=None,
    method="spades",
    memory=None,
    minlen=75,
    screen_accessions=None,
    screen_urls=None,
    mincontiglen=500,
    workdir=None,
    sourdb=None,
    mincovpct=5,
    debug=False,
    quiet=False,
    **kwargs,
):
    """Run the whole AAFTF pipeline, skipping steps whose output file already exists."""
    # passed to every step that has an option of the same name
    shared = {"cpus": cpus, "memory": memory, "workdir": workdir, "debug": debug, "quiet": quiet}
    trimmed = [basename + "_1P.fastq.gz", basename + "_2P.fastq.gz" if right else None]
    filtered = [basename + "_filtered_1.fastq.gz", basename + "_filtered_2.fastq.gz" if right else None]

    _run_step(trim, "trim", trimmed[0], shared, left=left, right=right, basename=basename, minlen=minlen)
    _run_step(aaftf_filter, "filter", filtered[0], shared, left=trimmed[0], right=trimmed[1], basename=basename, screen_accessions=screen_accessions, screen_urls=screen_urls)

    assembly = basename + f".{method}.fasta"
    _run_step(assemble, "assemble", assembly, shared, left=filtered[0], right=filtered[1], out=assembly, method=method, tmpdir=tmpdir, assembler_args=assembler_args)

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
        left=filtered[0],
        right=filtered[1],
        phylum=phylum,
        sourdb=sourdb,
        mincovpct=mincovpct,
    )

    rmdup_file = basename + ".rmdup.fasta"
    _run_step(rmdup, "rmdup", rmdup_file, shared, input=sourpurge_file, out=rmdup_file, minlen=mincontiglen)

    polish_file = basename + ".polish.fasta"
    _run_step(polish, "polish", polish_file, shared, infile=rmdup_file, outfile=polish_file, left=filtered[0], right=filtered[1])

    final_file = basename + ".final.fasta"
    _run_step(aaftf_sort, "sort", final_file, shared, input=polish_file, out=final_file, minlen=mincontiglen)

    assess.run(**_step_kwargs("assess", shared, input=final_file))


def _run_step(module, name, output, shared, **step_options):
    """Run ``module.run()`` for step ``name`` unless ``output`` already exists, then check that it does."""
    if check_file(output):
        logger.info(f"AAFTF {name} output found: {output}")
        return
    module.run(**_step_kwargs(name, shared, **step_options))
    if not check_file(output):
        raise RuntimeError(f"AAFTF {name} failed: {output} is missing or empty")


def _step_kwargs(name, shared, **step_options):
    """Return the run() keyword arguments for step ``name``: its CLI defaults, overridden by the pipeline's options.

    ``shared`` options are only passed to steps that have them, and a ``None``
    shared value (e.g. no pipeline ``--memory``) keeps the step's own default.
    Shared values take the type of the step's default (assemble's ``--memory`` is a string).
    """
    kwargs = dict(_subcommand_defaults()[name])
    for key, value in shared.items():
        if key in kwargs and value is not None:
            default = kwargs[key]
            kwargs[key] = type(default)(value) if type(default) in (int, str) else value
    kwargs.update(step_options)
    kwargs["pipe"] = True
    return kwargs


@functools.cache
def _subcommand_defaults():
    """Return ``{subcommand: {dest: default}}`` read from the AAFTF subcommand parsers."""
    from aaftf._menu import register_subcommands  # _menu imports this module, so import it when first needed

    subparsers = register_subcommands(ap.ArgumentParser())
    return {name: {action.dest: action.default for action in parser._actions if action.dest != "help"} for name, parser in subparsers.choices.items()}
