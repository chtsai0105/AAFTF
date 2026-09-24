"""Run NCBI routines to identify contaminant and vectior contigs.

The contaminants are presumably sequences that were not screened out
in the filter step.

This uses NCBI fcs tool for screening which relies on a singularity
engine installed

The default libraries for screening are located in resources.py
and include common Euk, Prok, and MITO contaminants.
"""

import logging
import shutil
import sys
from pathlib import Path

from aaftf.resources import FCSADAPTOR
from aaftf.utility import cleanup_workdir, make_workdir, require_databases, run_cmd

__all__ = ["run"]


logger = logging.getLogger(__name__)


def run(
    infile,
    outfile,
    container_engine="singularity",
    workdir=None,
    image=None,
    prok=False,
    fcs_script=None,
    debug=False,
    **kwargs,
):
    """Perform vector trimming via the fcs screening tool."""
    containerengine = container_engine
    infilename = Path(infile).resolve().name
    tax = "--euk"
    if prok:
        tax = "--prok"
    workdir, custom_workdir = make_workdir(workdir, "fcsscreen")

    # the wrapper script and singularity image come from `AAFTF download` unless given/found on PATH
    fcsexe = fcs_script or shutil.which("run_fcsadaptor.sh")
    needed = ([] if fcsexe else ["fcs_script"]) + (["fcs_image"] if containerengine == "singularity" and image is None else [])
    found = dict(zip(needed, require_databases(needed, hint="or pass --fcs_script PATH / --image PATH")))
    fcsexe = fcsexe or found["fcs_script"]

    if containerengine == "singularity":
        image = image or found["fcs_image"]
        if shutil.which("singularity") is None and shutil.which("apptainer") is None:
            logger.error("--container_engine singularity requires 'singularity' or 'apptainer' on PATH.")
            sys.exit(1)
    elif containerengine == "docker":
        # docker image reference (registry:tag), not a local file; docker itself
        # resolves/pulls it, so no download step is needed here.
        if image is None:
            image = FCSADAPTOR["DOCKERIMAGE"] % (FCSADAPTOR["VERSION"])
        if shutil.which("docker") is None:
            logger.error("--container_engine docker requires 'docker' on PATH.")
            sys.exit(1)

    cmd = [fcsexe, "--fasta-input", infile, "--output-dir", workdir, tax, "--container-engine", containerengine, "--image", image]
    run_cmd(cmd, debug)

    Path(workdir, "cleaned_sequences").mkdir()
    cleanresult = str(Path(workdir, "cleaned_sequences", infilename))
    if debug:
        logger.info(f"copy from: {cleanresult} -> {outfile}")
    Path(cleanresult).rename(outfile)
    fcsreport = str(Path(workdir, "fcs_adaptor_report.txt"))
    with open(fcsreport) as fh:
        logger.info("FCS report:")
        for line in fh:
            print(line, end="")
    # make a copy of the report to show
    Path(fcsreport).rename(outfile + ".fcs_adaptor_report.txt")
    # cleanup after running
    cleanup_workdir(workdir, debug, custom_workdir)
