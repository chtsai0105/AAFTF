#!/usr/bin/env python3
"""AAFTF main framework for submodule run."""

# note structure of code taken from poretools
# https://github.com/arq5x/poretools/blob/master/poretools/poretools_main.py

import argparse as ap
import logging
import sys

# AAFTF imports
from AAFTF import __version__
from AAFTF._menu import SUBCOMMAND_REGISTRARS
from AAFTF.utility import CustomHelpFormatter, available_cpus, getRAM, setup_logging

logger = logging.getLogger("AAFTF.main")

myversion = __version__


def _cap_to_available(args):
    """Lower -c/--cpus and -m/--memory to what this machine/job can provide, with a warning."""
    cpus = getattr(args, "cpus", None)
    if cpus and cpus > (avail_cpus := available_cpus()):
        logger.warning(f"-c/--cpus {cpus} is more than the {avail_cpus} CPUs available to this job; using {avail_cpus}")
        args.cpus = avail_cpus

    memory = getattr(args, "memory", None)
    if memory is not None and float(memory) > (avail_ram := getRAM()):
        capped = max(int(avail_ram), 1)
        logger.warning(f"-m/--memory {memory} GB is more than the {avail_ram:g} GB of RAM available; using {capped} GB")
        args.memory = type(memory)(capped)  # keep the option's type (int, or str for assemble)


def main():
    """Present the main AAFTF module submenus."""
    #########################################
    # create the top-level parser
    #########################################
    parser = ap.ArgumentParser(prog="AAFTF", formatter_class=CustomHelpFormatter)
    parser.add_argument("-V", "--version", help="Installed AAFTF version", action="version", version="%(prog)s " + str(myversion))

    subparsers = parser.add_subparsers(title="[sub-commands]", dest="command")

    #########################################
    # create the individual tool parsers
    #########################################
    for register in SUBCOMMAND_REGISTRARS:
        register(subparsers)

    # process the arguments now
    # if no args then print help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # each subcommand parser binds its own subtool's run() via
    # parser_x.set_defaults(func=<module>.run) in AAFTF/_menu.py; a bare
    # "AAFTF" invocation with unrecognized/no subcommand leaves func unset.
    if not hasattr(args, "func"):
        parser.print_help(sys.stderr)
        sys.exit(1)

    setup_logging(debug=getattr(args, "debug", False), quiet=getattr(args, "quiet", False))
    try:
        logger.info(f"Running AAFTF v{myversion}")
        _cap_to_available(args)
        args.func(**vars(args))
    except OSError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise


if __name__ == "__main__":
    main()
