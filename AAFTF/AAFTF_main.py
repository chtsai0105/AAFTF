#!/usr/bin/env python3
"""AAFTF main framework for submodule run."""

# note structure of code taken from poretools
# https://github.com/arq5x/poretools/blob/master/poretools/poretools_main.py

import argparse as ap
import sys

import AAFTF.assemble as assemble
import AAFTF.assess as assess
import AAFTF.check_dependencies as check_dependencies
import AAFTF.depth as depth
import AAFTF.download as download
import AAFTF.fcs_gx_purge as fcs_gx_purge
import AAFTF.fcs_screen as fcs_screen
import AAFTF.filter as aaftf_filter
import AAFTF.fix_tbl as fix_tbl
import AAFTF.mito as mito
import AAFTF.pipeline as pipeline
import AAFTF.polish as polish
import AAFTF.rmdup as rmdup
import AAFTF.sort as aaftf_sort
import AAFTF.sourpurge as sourpurge
import AAFTF.trim as trim
import AAFTF.vecscreen as vecscreen

# AAFTF imports
from AAFTF.__version__ import __version__, get_version  # noqa: E402
from AAFTF._menu import SUBCOMMAND_REGISTRARS
from AAFTF.utility import CustomHelpFormatter, status

try:
    myversion = get_version()
except Exception:
    # Fallback to hardcoded version if git isn't available
    myversion = __version__


ALIAS_MAP = {
    "trim_reads": "trim",
    "read_trim": "trim",
    "filter_reads": "filter",
    "read_filter": "filter",
    "asm": "assemble",
    "spades": "assemble",
    "vectorscreen": "vecscreen",
    "vector_blast": "vecscreen",
    "ncbi_fcs": "fcs_screen",
    "ncbi_fcs-screen": "fcs_screen",
    "ncbi_fcs-gx": "fcs_gx_purge",
    "ncbi_fcs_gx": "fcs_gx_purge",
    "gx": "fcs_gx_purge",
    "purge": "sourpurge",
    "dedup": "rmdup",
    "pilon": "polish",
    "polca": "polish",
    "stats": "assess",
    "fix": "fix_tbl",
    "mito_asm": "mito",
    "mitochondria": "mito",
    "coverage": "depth",
    "cov": "depth",
    "configure": "download",
    "install": "download",
    "download_db": "download",
    "setup": "download",
    "check_deps": "check_dependencies",
    "deps": "check_dependencies",
}

def run_subtool(parser, args):
    """Run the subtool of the AAFTF pipeline."""
    command = ALIAS_MAP.get(args.command, args.command)
    if command == "trim":
        submodule = trim
    elif command == "filter":
        submodule = aaftf_filter
    elif command == "assemble":
        submodule = assemble
    elif command == "vecscreen":
        submodule = vecscreen
    elif command == "fcs_screen":
        submodule = fcs_screen
    elif command == "fcs_gx_purge":
        submodule = fcs_gx_purge
    elif command == "sourpurge":
        submodule = sourpurge
    elif command == "rmdup":
        submodule = rmdup
    elif command == "polish":
        submodule = polish
    elif command == "assess":
        submodule = assess
    elif command == "sort":
        submodule = aaftf_sort
    elif command == "pipeline":
        submodule = pipeline
    elif command == "mito":
        submodule = mito
    elif command == "fix_tbl":
        submodule = fix_tbl
    elif command == "depth":
        submodule = depth
    elif command == "download":
        submodule = download
    elif command == "check_dependencies":
        submodule = check_dependencies
    else:
        parser.parse_args("")
        return
    # run the chosen submodule.
    submodule.run(parser, args)


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

    # done with menu options
    # set defaults
    parser.set_defaults(func=run_subtool)

    # process the arguments now
    # if no args then print help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    try:
        status(f"Running AAFTF v{myversion}")
        args.func(parser, args)
    except OSError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise


if __name__ == "__main__":
    main()
