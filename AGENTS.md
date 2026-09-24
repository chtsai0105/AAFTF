# AAFTF Development Guidelines

This document provides guidelines for agentic coding agents working on the AAFTF (Automatic Assembly For The Fungi) codebase.

## Project Overview

AAFTF is a Python-based bioinformatics toolkit for automated genome assembly, cleanup, mitochondrial genome assembly, and polishing focused on fungal genomes. The package uses external bioinformatics tools and provides a unified interface for genome assembly workflows.

## Build and Development Commands

### Installation and Setup
```bash
# Install in development mode
pip install -e .

# Install with dependencies from conda
conda create -n aaftf -c bioconda "python>=3.9" bbmap trimmomatic bowtie2 bwa sourmash \
    blast minimap2 spades megahit novoplasty biopython fastp pypolca polypolish nextpolish2 unicycler
pip install AAFTF
```

### Code Quality and Linting
```bash
# Run ruff for code style, linting, and import sorting (replaces flake8, isort, pyupgrade)
ruff check AAFTF/
ruff format AAFTF/

# Run ruff with automatic fixes
ruff check --fix AAFTF/

# Run pydocstyle for documentation standards
pydocstyle --convention=google AAFTF/

# Run codespell for spelling checks
codespell AAFTF/ --ignore-words-list=nd,reacher,thist,ths,ure,referenc,wile,nin,pid

# All pre-commit hooks (recommended)
pre-commit run --all-files
```

### Testing
```bash
# No formal pytest setup exists. Tests are shell scripts in tests/ directory
# Run individual test workflows:
bash tests/00_spades_1trim.sh
bash tests/01_mito.sh
bash tests/03_BUSCO.sh

# Basic import test
python3 -c "import AAFTF; print('Import successful')"
```

## Code Style Guidelines

### Import Organization
- Standard library imports first, then third-party, then local AAFTF imports
- Use explicit imports (avoid `from module import *`)
- Group related imports together
- Import AAFTF modules using relative imports within the package

```python
# Standard library
import os
import sys
import subprocess

# Third-party
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from packaging.version import Version

# Local AAFTF imports
from AAFTF.utility import status, checkfile, execute
from AAFTF.resources import DB_Links
```

### Function and Variable Naming
- Use `snake_case` for function and variable names
- Use `UPPER_SNAKE_CASE` for constants
- Function names should be descriptive verbs (`run_assembly`, `check_file`, `calculate_stats`)
- Variable names should be clear and concise (`input_file`, `output_dir`, `cpu_count`)

### Documentation and Docstrings
- Use Google-style docstrings
- Every module should have a module docstring explaining its purpose
- Every public function should have a docstring with description, arguments, and returns

```python
def calculate_n50(contig_lengths):
    """Calculate N50 statistic from a list of contig lengths.

    Args:
        contig_lengths: List of integers representing contig lengths

    Returns:
        Integer N50 value
    """
```

### Error Handling
- Use try/except blocks for external tool execution
- Provide informative error messages with context
- Report progress with a module-level `logger = logging.getLogger(__name__)` (`logger.info`, `logger.warning`, `logger.error`); `AAFTF_main.main()` configures it via `setup_logging()` from `-q/--quiet` and `-v/--debug`. Don't add "ERROR:"/"WARNING:" prefixes or leading/trailing whitespace to messages — the formatter adds the level name
- Handle file existence checks with `checkfile()` utility

### Code Structure
- Each AAFTF subcommand has its own module (trim.py, assemble.py, filter.py, etc.)
- Main entry point is through `AAFTF_main.py`
- Common utilities are in `utility.py`
- Resource URLs and constants are in `resources.py`
- Each module should have a `run(parser, args)` function for CLI integration

### CLI Argument Parsing
- Use `argparse` for command-line interface
- Follow existing patterns for common arguments (`-c/--cpus`, `-m/--memory`, `-o/--out`)
- Provide help text for all arguments
- Use consistent naming across subcommands

### External Tool Integration
- Check tool availability with `which_path()` utility
- Use `subprocess.run()` for external command execution
- Handle tool version differences with `packaging.version.Version`
- Provide fallback options for different installation methods (conda, homebrew, etc.)

### File I/O Patterns
- Support both compressed (.gz) and uncompressed files
- Use `checkfile()` to validate input files
- Use `safe_remove()` for file cleanup. Create a subcommand's working directory with `workdir, custom_workdir = make_workdir(workdir, "<name>")` and remove it with `cleanup_workdir(workdir, debug, custom_workdir)` (it never deletes a user-supplied `--workdir`)
- Other shared helpers in `AAFTF/utility.py` — use these instead of re-implementing: `download_file()` (never leaves a partial file behind), `aaftf_db_dir()`, `require_tools()`, `run_cmd()` / `execute()`, `align_to_sorted_bam()`, `concat_files()`, `open_maybe_gz()`, `next_step_name()`, `basename_from_reads()`, `filter_fasta()` / `write_fasta()`, `calc_nx()`
- Handle file paths with `os.path` operations for cross-platform compatibility

### Logging and Output
- Use `logger.info/warning/error` for user-facing messages
- Use `printCMD()` to show commands being executed
- Provide progress feedback for long-running operations
- Use debug flags for verbose output during development

### Performance Considerations
- Respect CPU and memory limits specified by user
- Use streaming operations for large files when possible
- Implement efficient file parsing with BioPython iterators
- Consider external tool memory requirements in parameter defaults

### Version Handling
- `AAFTF/__init__.py` resolves `AAFTF.__version__` via `importlib.metadata.version("AAFTF")`, falling back to `"0.0.0+unknown"` if the package isn't installed
- The version string itself is set at build time by `hatch-vcs` from git tags (see `[tool.hatch.version]` in `pyproject.toml`)
- Main application (`AAFTF_main.py`) imports it with `from AAFTF import __version__`
- Version is displayed both via `--version` flag and at application startup
- Maintains PEP 440 compatibility for Python packaging standards

## Testing Guidelines

Since AAFTF primarily integrates external bioinformatics tools, testing focuses on:
- Integration workflows with real data (shell scripts in tests/)
- Import and basic functionality testing
- CLI argument parsing validation
- File format compatibility (FASTA/FASTQ parsing)

## Common Patterns

### Subcommand Module Structure
```python
"""Module description."""

import logging
import sys

from AAFTF.utility import checkfile, run_cmd

logger = logging.getLogger(__name__)


def run(input, debug=False, pipe=False, **kwargs):
    """Main entry point for subcommand."""
    # Validate inputs
    if not checkfile(input):
        logger.error(f"Input file not found: {input}")
        sys.exit(1)

    logger.info("Starting processing...")

    # Execute external tools (prints the command; stderr shown only with --debug)
    run_cmd(["tool", "--input", input], debug)

    logger.info("Processing complete")
```

This project follows bioinformatics best practices with emphasis on reproducibility, proper tool integration, and clear user feedback.
