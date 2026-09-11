"""Run the fcs-gx tool to look for contaminants."""

import shutil
import subprocess
import sys
import uuid
from pathlib import Path

from Bio import SeqIO

from AAFTF.utility import SafeRemove, checkfile, fastastats, printCMD, status

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
    if not workdir:
        workdir = f"aaftf-fcsgx_{str(uuid.uuid4())[:8]}"
    if not Path(workdir).exists():
        Path(workdir).mkdir()

    # parse database locations
    if not db or not Path(f"{db}.gxi").is_file():
        status(f"{db}.gxi not found, needs to have setup the fcs_gx db - https://github.com/ncbi/fcs/wiki/FCS-GX")
        sys.exit(1)

    numSeqs, assemblySize = fastastats(str(Path(input)))
    status(f"Assembly is {numSeqs:,} contigs and {assemblySize:,} bp")

    # now filter for taxonomy with sourmash lca classify
    status("Running fcs_gx to get taxonomy classification for each contig")

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
        status(f"fcs_gx did not produce file {fcsgxTSV}")
        return
    with open(fcsgxTSV) as fcsgx_out:
        for line in fcsgx_out:
            if not line.startswith("#"):
                line = line.strip()
                cols = line.split("\t")
                Seq2Drop[cols[0]] = 1

    # drop contigs from taxonomy before calculating coverage
    status(f"Dropping {len(Seq2Drop)} contigs from fcs-gx taxonomy screen")
    with open(outfile, "w") as ofh:
        for record in SeqIO.parse(input, "fasta"):
            if record.id not in Seq2Drop:
                SeqIO.write(record, ofh, "fasta")

    if debug:
        print("Contigs dropped due to taxonomy: {:}".format(",".join(Seq2Drop)))

    numSeqs, assemblySize = fastastats(outfile)
    status(f"fcs-gx assembly is {numSeqs:,} contigs and {assemblySize:,} bp")
    if "_" in outfile:
        nextOut = outfile.split("_")[0] + ".rmdup.fasta"
    elif "." in outfile:
        nextOut = outfile.split(".")[0] + ".rmdup.fasta"
    else:
        nextOut = f"{outfile}.rmdup.fasta"

    if checkfile(fcsgxTSV):
        outbase = Path(outfile).name
        basedir = str(Path(outfile).resolve().parent)
        if "." in outbase:
            outbase = outbase.rsplit(".", 1)[0]
        shutil.copy(fcsgxTSV, str(Path(basedir, f"{outbase}.fcs_gx-taxonomy.tsv")))

    if not debug:
        SafeRemove(workdir)

    if not pipe:
        status(f"Your next command might be:\n\tAAFTF rmdup -i {outfile} -o {nextOut}\n")
