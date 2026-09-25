"""Check that AAFTF's external tool and Python package dependencies are installed."""

import importlib
import logging
import shutil
from typing import Any

__all__ = ["REQUIRED_TOOLS", "OPTIONAL_TOOLS", "REQUIRED_PYTHON_PACKAGES", "run"]


# external binaries used via subprocess across the AAFTF modules
REQUIRED_TOOLS = {
    "samtools": "filter, sourpurge, rmdup, polish, depth",
    "bwa": "filter, sourpurge, polish, depth",
    "minimap2": "mito, rmdup, depth, polish (nextpolish2, racon)",
    "racon": "polish (racon)",
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
    "reformat.sh": "trim, filter, mito (--subsample)",
    "shuffle.sh": "trim, filter",
    "pypolca": "polish (pypolca)",
    "freebayes": "polish (pypolca)",
    "polypolish": "polish (polypolish)",
    "nextPolish2": "polish (nextpolish2)",
    "yak": "polish (nextpolish2)",
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

logger = logging.getLogger(__name__)


def run(**kwargs: Any) -> None:
    """Check whether AAFTF's external tool and Python package dependencies are installed.

    Prints an OK/MISSING line for every required and optional tool and required package.

    Args:
        **kwargs: Parsed CLI attributes (``command``, ``func``, ``debug``, ...); ignored.

    Raises:
        FileNotFoundError: If any required tool or Python package is missing.
    """
    missing_required = _print_tool_table("Checking required external tools...", _check_tools(REQUIRED_TOOLS))
    print()
    missing_optional = _print_tool_table("Checking optional external tools...", _check_tools(OPTIONAL_TOOLS))

    print()
    logger.info("Checking required Python packages...")
    missing_packages = []
    for dist_name, found in _check_python_packages(REQUIRED_PYTHON_PACKAGES):
        if found:
            print(f"  [OK]      {dist_name}")
        else:
            print(f"  [MISSING] {dist_name}")
            missing_packages.append(dist_name)

    print()
    if missing_required or missing_packages:
        if missing_optional:
            logger.info(f"NOTE: {len(missing_optional)} optional tool(s) also missing (only needed for specific subcommands).")
        raise FileNotFoundError(f"{len(missing_required)} required tool(s) and {len(missing_packages)} required package(s) missing.")
    else:
        logger.info("All required dependencies are installed.")
        if missing_optional:
            logger.info(f"NOTE: {len(missing_optional)} optional tool(s) missing (only needed for specific subcommands).")


def _print_tool_table(title: str, results: list[tuple[str, str | None, str]]) -> list[str]:
    """Log ``title`` and print an OK/MISSING line for each checked tool.

    Args:
        title: Heading logged before the table.
        results: ``(tool, path, used_by)`` tuples from ``_check_tools``; ``path`` is None if not found.

    Returns:
        Names of the tools that were not found.
    """
    logger.info(title)
    missing = []
    for tool, path, used_by in results:
        if path:
            print(f"  [OK]      {tool:<20} {path}")
        else:
            print(f"  [MISSING] {tool:<20} required by: {used_by}")
            missing.append(tool)
    return missing


def _check_tools(tools: dict[str, str]) -> list[tuple[str, str | None, str]]:
    """Look up each tool on ``PATH``.

    Args:
        tools: Mapping of executable name to the subcommands that use it.

    Returns:
        ``(tool, path, used_by)`` tuples sorted by tool name; ``path`` is None if not on ``PATH``.
    """
    results = []
    for tool, used_by in sorted(tools.items()):
        path = shutil.which(tool)
        results.append((tool, path, used_by))
    return results


def _check_python_packages(packages: dict[str, str]) -> list[tuple[str, bool]]:
    """Try importing each Python package.

    Args:
        packages: Mapping of import name to distribution name.

    Returns:
        ``(distribution name, importable)`` tuples sorted by import name.
    """
    results = []
    for module_name, dist_name in sorted(packages.items()):
        try:
            importlib.import_module(module_name)
            results.append((dist_name, True))
        except ImportError:
            results.append((dist_name, False))
    return results
