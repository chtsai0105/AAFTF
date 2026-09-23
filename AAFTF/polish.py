"""Runs short/long-read polishing of an assembly.

Three polishing engines are supported via --method:
  - pypolca:    single-pass POLCA-style polishing (Illumina short reads)
  - polypolish: alignment-filtering short-read polisher (Illumina short reads)
  - nextpolish2: repeat-aware polishing of HiFi assemblies using a short-read
                 k-mer (yak) database (requires --longreads HiFi + short reads)
"""

import shutil
import subprocess
import sys
import uuid
from pathlib import Path

from AAFTF.utility import align_to_sorted_bam, cleanup_workdir, printCMD, run_cmd, status


def run(
    infile,
    outfile=None,
    method="polypolish",
    memory=16,
    cpus=1,
    left=None,
    right=None,
    longreads=None,
    workdir=None,
    debug=False,
    pipe=False,
    **kwargs,
):
    """Execute polishing step provided with reads and a contig assembly FASTA file."""
    method = method.lower()
    status(f"calling {method} with input memory {memory}GB and num cpus {cpus}")

    forReads, revReads = (None,) * 2
    if left:
        forReads = str(Path(left).resolve())
    if right:
        revReads = str(Path(right).resolve())
    if longreads:
        longreads = str(Path(longreads).resolve())

    if method == "racon" and not longreads:
        status("Unable to locate long read FASTQ raw reads, pass via -lr or --longreads")
        sys.exit(1)
    if method == "nextpolish2" and not longreads:
        status("Unable to locate long read FASTQ raw reads, pass via -lr or --longreads (nextpolish2 requires HiFi long reads)")
        sys.exit(1)
    if method in ("pypolca", "masurca", "polypolish", "nextpolish2") and not forReads:
        status("Unable to locate FASTQ raw reads, pass via -l,--left and/or -r,--right")
        sys.exit(1)

    custom_workdir = bool(workdir)
    if not workdir:
        workdir = f"aaftf-polish_{str(uuid.uuid4())[:8]}"
    if not Path(workdir).exists():
        Path(workdir).mkdir()

    # Output file
    polishedFasta = outfile
    if not polishedFasta:
        fbasename = Path(infile).name.split(".f")[0]
        polishedFasta = f"{fbasename}.polished.fasta"

    polish_log = f"{method}.log"

    # Preflight: make sure the external tools this method needs are on PATH.
    required_exes = {
        "polypolish": ["bwa", "polypolish"],
        "pypolca": ["bwa", "samtools", "freebayes", "pypolca"],
        "nextpolish2": ["minimap2", "samtools", "yak", "nextPolish2"],
        "racon": ["minimap2", "racon"],
    }.get(method, [])
    missing = [exe for exe in required_exes if shutil.which(exe) is None]
    if missing:
        status(f"ERROR: required executable(s) not found on PATH for --method {method}: {', '.join(missing)}")
        status("Install the missing tool(s) (e.g. `pixi add polypolish` / `conda install -c bioconda polypolish`) and ensure the correct environment is activated.")
        sys.exit(1)

    if method == "polypolish":
        ret, out_path = run_polypolish(infile, forReads, revReads, cpus, workdir, polish_log, debug)
    elif method == "pypolca":
        ret, out_path = run_pypolca(infile, forReads, revReads, cpus, memory, workdir, polish_log, polishedFasta)
    elif method == "nextpolish2":
        ret, out_path = run_nextpolish2(infile, forReads, revReads, longreads, cpus, workdir, polish_log, debug)
    elif method == "racon":
        ret, out_path = run_racon(infile, longreads, cpus, workdir, polish_log, debug)
    else:
        status(f"Unknown polishing method: {method}")
        sys.exit(1)

    # Validate the polisher's output and copy it to the requested destination
    # — done once here rather than duplicated in every run_<method>() function.
    if ret != 0 or not Path(out_path).exists() or Path(out_path).stat().st_size == 0:
        status(f"ERROR: {method} failed (exit {ret}); check log: {Path(workdir, polish_log)}")
        sys.exit(1)
    shutil.copyfile(out_path, polishedFasta)
    status("AAFTF polish completed.")
    status(f"{method} polished assembly: {polishedFasta}")

    nextOut = _derive_next_out(polishedFasta)

    cleanup_workdir(workdir, debug, custom_workdir)

    if not pipe:
        status("Your next command might be:\n" + f"\tAAFTF sort -i {polishedFasta} -o {nextOut}\n")


def _derive_next_out(polishedFasta):
    """Derive the suggested next-step (sort) output filename from polishedFasta."""
    if "_" in polishedFasta:
        return polishedFasta.split("_")[0] + ".final.fasta"
    elif "." in polishedFasta:
        return polishedFasta.split(".")[0] + ".final.fasta"
    return polishedFasta + ".final.fasta"


def run_polypolish(infile, forReads, revReads, cpus, workdir, polish_log, debug):
    """Polish with Polypolish.

    Workflow: index the assembly, align each read file separately with
    ``bwa mem -a`` (reporting all alignments, required by Polypolish), filter
    the resulting SAM pairs by insert size, then polish.

    Returns (returncode, out_path). Validation, copying to the final
    destination, and status reporting are all handled centrally by run().
    """
    if not revReads:
        status("ERROR: --method polypolish requires paired reads (-l/--left and -r/--right)")
        sys.exit(1)

    assembly = str(Path(workdir, Path(infile).name))
    shutil.copyfile(infile, assembly)
    asm_name = Path(assembly).name

    bwa_index = ["bwa", "index", asm_name]
    run_cmd(bwa_index, debug, cwd=workdir)

    sam1, sam2 = "alignments_1.sam", "alignments_2.sam"
    for reads, sam_out in ((forReads, sam1), (revReads, sam2)):
        bwa_cmd = ["bwa", "mem", "-t", str(cpus), "-a", asm_name, reads]
        with open(str(Path(workdir, sam_out)), "w") as out_fh:
            run_cmd(bwa_cmd, debug, cwd=workdir, stdout=out_fh)

    filt1, filt2 = "filtered_1.sam", "filtered_2.sam"
    filter_cmd = ["polypolish", "filter", "--in1", sam1, "--in2", sam2, "--out1", filt1, "--out2", filt2]
    run_cmd(filter_cmd, debug, cwd=workdir)

    polish_cmd = ["polypolish", "polish", asm_name, filt1, filt2]
    printCMD(polish_cmd)
    out_path = str(Path(workdir, "polypolish_corrected.fasta"))
    with open(str(Path(workdir, polish_log)), "w") as logfile, open(out_path, "w") as out_fh:
        ret = subprocess.run(polish_cmd, cwd=workdir, stdout=out_fh, stderr=logfile)
    return ret.returncode, out_path


def run_pypolca(infile, forReads, revReads, cpus, memory, workdir, polish_log, polishedFasta):
    """Polish with pypolca (POLCA algorithm reimplemented in Python; runs bwa+samtools+freebayes internally).

    Copies pypolca's ``.vcf``/``.report`` sidecar outputs alongside
    ``polishedFasta`` itself (as ``{polishedFasta}.vcf`` /
    ``{polishedFasta}.pypolca_report.txt``), since those are specific to this
    method. Returns (returncode, out_path) for the polished FASTA; validating
    and copying that primary output is handled centrally by run().
    """
    memperthread = int(memory / cpus)
    if memperthread == 0:
        # minimum should be 1Gb
        memperthread = "1"
    memperthread = f"{memperthread}G"

    pypolca_prefix = "pypolca"
    pypolca_outdir = str(Path(workdir, "pypolca_out"))
    pypolca_cmd = ["pypolca", "run", "-a", str(Path(infile).resolve()), "-1", forReads]
    if revReads:
        pypolca_cmd.extend(["-2", revReads])
    pypolca_cmd.extend(["-t", str(cpus), "-o", pypolca_outdir, "-p", pypolca_prefix, "-m", memperthread, "-f"])
    printCMD(pypolca_cmd)
    with open(str(Path(workdir, polish_log)), "w") as logfile:
        ret = subprocess.run(pypolca_cmd, cwd=workdir, stderr=logfile, stdout=logfile)
    out_path = str(Path(pypolca_outdir, f"{pypolca_prefix}_corrected.fasta"))
    if ret.returncode == 0:
        vcf_src = str(Path(pypolca_outdir, f"{pypolca_prefix}.vcf"))
        report_src = str(Path(pypolca_outdir, f"{pypolca_prefix}.report"))
        if Path(vcf_src).exists():
            shutil.copyfile(vcf_src, f"{polishedFasta}.vcf")
        if Path(report_src).exists():
            shutil.copyfile(report_src, f"{polishedFasta}.pypolca_report.txt")
    return ret.returncode, out_path


def run_nextpolish2(infile, forReads, revReads, longreads, cpus, workdir, polish_log, debug):
    """Polish with NextPolish2.

    Workflow: build a short-read k-mer (yak) database, map HiFi long reads to
    the assembly with minimap2, then run nextPolish2's repeat-aware
    correction using the HiFi alignments plus the yak k-mer database.

    Returns (returncode, out_path). Validation, copying to the final
    destination, and status reporting are all handled centrally by run().
    """
    assembly = str(Path(workdir, Path(infile).name))
    shutil.copyfile(infile, assembly)
    asm_name = Path(assembly).name

    yak_db = "shortreads.yak"
    yak_cmd = ["yak", "count", "-k31", "-b37", "-t", str(cpus), "-o", yak_db, forReads]
    if revReads:
        yak_cmd.append(revReads)
    run_cmd(yak_cmd, debug, cwd=workdir)

    hifi_bam = "hifi.map.bam"
    minimap_cmd = ["minimap2", "-ax", "map-hifi", "-t", str(cpus), asm_name, longreads]
    align_to_sorted_bam(minimap_cmd, str(Path(workdir, hifi_bam)), cpus, cwd=workdir, debug=debug)

    out_fasta = "nextpolish2_corrected.fasta"
    nextpolish2_cmd = ["nextPolish2", "-t", str(cpus), "-o", out_fasta, hifi_bam, asm_name, yak_db]
    printCMD(nextpolish2_cmd)
    out_path = str(Path(workdir, out_fasta))
    with open(str(Path(workdir, polish_log)), "w") as logfile:
        ret = subprocess.run(nextpolish2_cmd, cwd=workdir, stderr=logfile, stdout=logfile)
    return ret.returncode, out_path


def run_racon(infile, longreads, cpus, workdir, polish_log, debug):
    """Polish with Racon.

    Workflow: map long reads to the assembly with minimap2 (producing PAF
    overlaps), then run racon using those overlaps to correct the assembly.

    Returns (returncode, out_path). Validation, copying to the final
    destination, and status reporting are all handled centrally by run().
    """
    assembly = str(Path(workdir, Path(infile).name))
    shutil.copyfile(infile, assembly)
    asm_name = Path(assembly).name

    overlaps = "overlaps.paf"
    minimap_cmd = ["minimap2", "-x", "map-ont", "-t", str(cpus), asm_name, longreads]
    with open(str(Path(workdir, overlaps)), "w") as out_fh:
        run_cmd(minimap_cmd, debug, cwd=workdir, stdout=out_fh)

    racon_cmd = ["racon", "-t", str(cpus), longreads, overlaps, asm_name]
    printCMD(racon_cmd)
    out_path = str(Path(workdir, "racon_corrected.fasta"))
    with open(str(Path(workdir, polish_log)), "w") as logfile, open(out_path, "w") as out_fh:
        ret = subprocess.run(racon_cmd, cwd=workdir, stdout=out_fh, stderr=logfile)
    return ret.returncode, out_path
