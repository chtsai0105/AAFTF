"""Run NCBI routines to identify contaminant and vectior contigs.

The contaminants are presumably sequences that were not screened out
in the filter step.

This uses NCBI fcs tool for screening which relies on a singularity
engine installed

The default libraries for screening are located in resources.py
and include common Euk, Prok, and MITO contaminants.
"""

import shutil
import sys
from pathlib import Path

from AAFTF.resources import FCSADAPTOR
from AAFTF.utility import aaftf_db_dir, cleanup_workdir, download_file, make_workdir, run_cmd, status


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
    DB = aaftf_db_dir(required=True)

    containerengine = container_engine
    infilename = Path(infile).resolve().name
    tax = "--euk"
    if prok:
        tax = "--prok"
    workdir, custom_workdir = make_workdir(workdir, "fcsscreen")

    fcsexe = fcs_script
    if fcsexe is None:
        fcsexe = shutil.which("run_fcsadaptor.sh")
    if fcsexe is None:
        fcsexe = str(Path(DB, "run_fcsadaptor.sh"))
        #  This will help download the fcs-adaptor shell script rather than re-implementing it here
        if not Path(fcsexe).exists():
            url = FCSADAPTOR["EXEURL"] % (FCSADAPTOR["VERSION"])
            download_file(url, fcsexe)
            Path(fcsexe).chmod(0o555)

    if containerengine == "singularity":
        # local SIF file: download once and cache under AAFTF_DB
        if image is None:
            image = str(Path(DB, FCSADAPTOR["SIFLOCAL"] % (FCSADAPTOR["VERSION"])))
            if not Path(image).exists():
                # SIFURL is a URL prefix, not a filesystem path — do not use pathlib here.
                url = "/".join([FCSADAPTOR["SIFURL"].rstrip("/"), FCSADAPTOR["VERSION"], FCSADAPTOR["SIF"]])
                download_file(url, image)
        if shutil.which("singularity") is None and shutil.which("apptainer") is None:
            status("ERROR: --container_engine singularity requires 'singularity' or 'apptainer' on PATH.")
            sys.exit(1)
    elif containerengine == "docker":
        # docker image reference (registry:tag), not a local file; docker itself
        # resolves/pulls it, so no download step is needed here.
        if image is None:
            image = FCSADAPTOR["DOCKERIMAGE"] % (FCSADAPTOR["VERSION"])
        if shutil.which("docker") is None:
            status("ERROR: --container_engine docker requires 'docker' on PATH.")
            sys.exit(1)

    cmd = [fcsexe, "--fasta-input", infile, "--output-dir", workdir, tax, "--container-engine", containerengine, "--image", image]
    run_cmd(cmd, debug)

    Path(workdir, "cleaned_sequences").mkdir()
    cleanresult = str(Path(workdir, "cleaned_sequences", infilename))
    if debug:
        status(f"copy from: {cleanresult} -> {outfile}")
    Path(cleanresult).rename(outfile)
    fcsreport = str(Path(workdir, "fcs_adaptor_report.txt"))
    with open(fcsreport) as fh:
        status("FCS report:")
        for line in fh:
            print(line, end="")
    # make a copy of the report to show
    Path(fcsreport).rename(outfile + ".fcs_adaptor_report.txt")
    # cleanup after running
    cleanup_workdir(workdir, debug, custom_workdir)
