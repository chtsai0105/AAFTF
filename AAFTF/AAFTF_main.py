#!/usr/bin/env python3
"""AAFTF main framework for submodule run."""

# note structure of code taken from poretools
# https://github.com/arq5x/poretools/blob/master/poretools/poretools_main.py

import argparse as ap
import sys

# AAFTF imports
from AAFTF._menu import SUBCOMMAND_REGISTRARS
from AAFTF._version import __version__
from AAFTF.utility import CustomHelpFormatter, status

myversion = __version__


def main():
    """Present the main AAFTF module submenus."""
    #########################################
    # create the top-level parser
    #########################################
    parser = ap.ArgumentParser(prog="AAFTF", formatter_class=CustomHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed AAFTF version", action="version", version="%(prog)s " + str(myversion))

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

    try:
        status(f"Running AAFTF v{myversion}")
        args.func(**vars(args))
    except OSError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise


if __name__ == "__main__":
    main()
