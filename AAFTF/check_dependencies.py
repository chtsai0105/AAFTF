"""Check that AAFTF's external tool and Python package dependencies are installed."""

import importlib
import sys

from AAFTF.utility import status, which_path

# external binaries used via subprocess across the AAFTF modules
REQUIRED_TOOLS = {
    "samtools": "filter, sourpurge, rmdup, polish, depth",
    "bwa": "filter, sourpurge, polish, depth",
    "minimap2": "mito, rmdup, depth",
    "bowtie2": "filter",
    "bowtie2-build": "filter",
    "blastn": "vecscreen",
    "makeblastdb": "vecscreen",
    "sourmash": "sourpurge",
    "mosdepth": "depth",
    "spades.py": "assemble",
    "megahit": "assemble",
    "unicycler": "assemble",
    "fastp": "trim",
    "trimmomatic": "trim",
    "java": "trim (trimmomatic)",
    "bbduk.sh": "trim, filter",
    "reformat.sh": "trim, filter",
    "shuffle.sh": "trim, filter",
    "pilon": "polish",
    "nextPolish": "polish",
    "polca.sh": "polish (masurca)",
    "NOVOPlasty.pl": "mito",
}

# optional tools only needed for specific, opt-in subcommands/container backends
OPTIONAL_TOOLS = {
    "singularity": "fcs_screen (--container_engine singularity)",
    "apptainer": "fcs_screen (--container_engine singularity)",
    "docker": "fcs_screen (--container_engine docker)",
    "run_fcsadaptor.sh": "fcs_screen",
    "run_gx.py": "fcs_gx_purge",
}

REQUIRED_PYTHON_PACKAGES = {
    "Bio": "biopython",
}


def _check_tools(tools):
    results = []
    for tool, used_by in sorted(tools.items()):
        path = which_path(tool)
        results.append((tool, path, used_by))
    return results


def _check_python_packages(packages):
    results = []
    for module_name, dist_name in sorted(packages.items()):
        try:
            importlib.import_module(module_name)
            results.append((dist_name, True))
        except ImportError:
            results.append((dist_name, False))
    return results


def _print_tool_table(title, results):
    status(title)
    missing = []
    for tool, path, used_by in results:
        if path:
            print(f"  [OK]      {tool:<20} {path}")
        else:
            print(f"  [MISSING] {tool:<20} required by: {used_by}")
            missing.append(tool)
    return missing


def run(parser, args):
    """Check whether AAFTF's external tool and Python package dependencies are installed."""
    missing_required = _print_tool_table("Checking required external tools...", _check_tools(REQUIRED_TOOLS))
    print()
    missing_optional = _print_tool_table("Checking optional external tools...", _check_tools(OPTIONAL_TOOLS))

    print()
    status("Checking required Python packages...")
    missing_packages = []
    for dist_name, found in _check_python_packages(REQUIRED_PYTHON_PACKAGES):
        if found:
            print(f"  [OK]      {dist_name}")
        else:
            print(f"  [MISSING] {dist_name}")
            missing_packages.append(dist_name)

    print()
    if missing_required or missing_packages:
        status(f"ERROR: {len(missing_required)} required tool(s) and {len(missing_packages)} required package(s) missing.")
        if missing_optional:
            status(f"NOTE: {len(missing_optional)} optional tool(s) also missing (only needed for specific subcommands).")
        sys.exit(1)
    else:
        status("All required dependencies are installed.")
        if missing_optional:
            status(f"NOTE: {len(missing_optional)} optional tool(s) missing (only needed for specific subcommands).")
