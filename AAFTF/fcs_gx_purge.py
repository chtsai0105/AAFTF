"""Run the fcs-gx tool to look for contaminants."""

import logging
import shutil
import subprocess
import sys
from pathlib import Path

from AAFTF.utility import checkfile, cleanup_workdir, fastastats, filter_fasta, make_workdir, next_step_name, printCMD

logger = logging.getLogger(__name__)

# logging - we may need to think about whether this has
# separate name for the different runfolder


# flake8: noqa: C901
def run(
    input,
    outfile,
    workdir=None,
    taxid=4890,
    db="/my_tmpfs/gxdb/all",
    debug=False,
    pipe=False,
    **kwargs,
):
    """Run NCBI fcs_gx routines to detect and remove contaminant contigs."""
    workdir, custom_workdir = make_workdir(workdir, "fcsgx")

    # parse database locations
    if not db or not Path(f"{db}.gxi").is_file():
        logger.info(f"{db}.gxi not found, needs to have setup the fcs_gx db - https://github.com/ncbi/fcs/wiki/FCS-GX")
        sys.exit(1)

    numSeqs, assemblySize = fastastats(str(Path(input)))
    logger.info(f"Assembly is {numSeqs:,} contigs and {assemblySize:,} bp")

    # now filter for taxonomy with sourmash lca classify
    logger.info("Running fcs_gx to get taxonomy classification for each contig")

    # python scripts/run_gx.py --bin-dir dist --gx-db /sw/db/gxdb --tax-id 4842 --fasta

    fcsgx_compute = [
        "run_gx.py",
        "--fasta",
        input,
        "--tax-id",
        f"{taxid}",  # taxid is a numeric
        "--gx-db",
        db,
        "--out-dir",
        workdir,
    ]
    printCMD(fcsgx_compute)
    fcs_log = "fcs_gx.log"
    with open(str(Path(workdir, fcs_log)), "w") as logfile:
        subprocess.run(fcsgx_compute, stderr=logfile)

    fname = Path(input).stem
    # output tsv:
    # seq_id	start_pos	end_pos	seq_len	action	div	agg_cont_cov	top_tax_name

    Seq2Drop = {}

    fcsgxTSV = str(Path(workdir, f"{fname}.{taxid}.fcs_gx_report.txt"))
    if not Path(fcsgxTSV).is_file():
        logger.info(f"fcs_gx did not produce file {fcsgxTSV}")
        return
    with open(fcsgxTSV) as fcsgx_out:
        for line in fcsgx_out:
            if not line.startswith("#"):
                line = line.strip()
                cols = line.split("\t")
                Seq2Drop[cols[0]] = 1

    # drop contigs from taxonomy before calculating coverage
    logger.info(f"Dropping {len(Seq2Drop)} contigs from fcs-gx taxonomy screen")
    numSeqs, assemblySize = filter_fasta(input, outfile, lambda seq_id: seq_id not in Seq2Drop)

    if debug:
        print("Contigs dropped due to taxonomy: {:}".format(",".join(Seq2Drop)))

    logger.info(f"fcs-gx assembly is {numSeqs:,} contigs and {assemblySize:,} bp")
    nextOut = next_step_name(outfile, ".rmdup.fasta")

    if checkfile(fcsgxTSV):
        outbase = Path(outfile).name
        basedir = str(Path(outfile).resolve().parent)
        if "." in outbase:
            outbase = outbase.rsplit(".", 1)[0]
        shutil.copy(fcsgxTSV, str(Path(basedir, f"{outbase}.fcs_gx-taxonomy.tsv")))

    cleanup_workdir(workdir, debug, custom_workdir)

    if not pipe:
        logger.info(f"Your next command might be:\nAAFTF rmdup -i {outfile} -o {nextOut}")
