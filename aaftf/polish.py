"""Polish an assembly with short and/or long reads.

Four polishing engines are supported via --method:
  - pypolca:    single-pass POLCA-style polishing (Illumina short reads)
  - polypolish: alignment-filtering short-read polisher (Illumina short reads)
  - nextpolish2: repeat-aware polishing of HiFi assemblies using a short-read
                 k-mer (yak) database (requires --longreads HiFi + short reads)
  - racon:      long-read consensus polishing from minimap2 overlaps (requires --longreads)
"""

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any

from aaftf.utility import align_to_sorted_bam, cleanup_workdir, make_workdir, next_step_name, print_cmd, require_tools, run_cmd

__all__ = ["run", "run_polypolish", "run_pypolca", "run_nextpolish2", "run_racon"]


logger = logging.getLogger(__name__)


def run(
    infile: str,
    outfile: str | None = None,
    method: str = "polypolish",
    memory: int = 16,
    cpus: int = 1,
    left: str | None = None,
    right: str | None = None,
    longreads: str | None = None,
    workdir: str | None = None,
    debug: bool = False,
    pipe: bool = False,
    **kwargs: Any,
) -> None:
    """Run the ``polish`` subcommand: polish an assembly FASTA with the chosen method.

    Checks the required reads and tools, runs the method's ``run_<method>()`` helper, validates
    its output and copies it to ``outfile``.

    Args:
        infile: Input assembly FASTA.
        outfile: Output polished FASTA; ``<input prefix>.polished.fasta`` if None.
        method: ``"polypolish"``, ``"pypolca"``, ``"nextpolish2"`` or ``"racon"`` (case-insensitive).
        memory: Total memory in GB (pypolca only).
        cpus: Number of threads.
        left: Left/forward short reads, or None.
        right: Right/reverse short reads, or None.
        longreads: Long-read FASTQ, or None.
        workdir: Working directory; a temporary one is created if None.
        debug: Keep the working directory and show command output when True.
        pipe: Suppress the "next command" hint when True.
        **kwargs: Other parsed CLI attributes (``command``, ``func``, ``quiet``); ignored.

    Raises:
        ValueError: If reads required by ``method`` are missing or ``method`` is unknown.
        RuntimeError: If the polisher fails or produces an empty output.
    """
    method = method.lower()
    logger.info(f"calling {method} with input memory {memory}GB and num cpus {cpus}")

    forward_reads, reverse_reads = (None,) * 2
    if left:
        forward_reads = str(Path(left).resolve())
    if right:
        reverse_reads = str(Path(right).resolve())
    if longreads:
        longreads = str(Path(longreads).resolve())

    if method == "racon" and not longreads:
        raise ValueError("Unable to locate long read FASTQ raw reads, pass via -lr or --longreads")
    if method == "nextpolish2" and not longreads:
        raise ValueError("Unable to locate long read FASTQ raw reads, pass via -lr or --longreads (nextpolish2 requires HiFi long reads)")
    if method in ("pypolca", "masurca", "polypolish", "nextpolish2") and not forward_reads:
        raise ValueError("Unable to locate FASTQ raw reads, pass via -l,--left and/or -r,--right")

    workdir, custom_workdir = make_workdir(workdir, "polish")

    # Output file
    polished_fasta = outfile
    if not polished_fasta:
        fbasename = Path(infile).name.split(".f")[0]
        polished_fasta = f"{fbasename}.polished.fasta"

    polish_log = f"{method}.log"

    # Preflight: make sure the external tools this method needs are on PATH.
    required_exes = {
        "polypolish": ["bwa", "polypolish"],
        "pypolca": ["bwa", "samtools", "freebayes", "pypolca"],
        "nextpolish2": ["minimap2", "samtools", "yak", "nextPolish2"],
        "racon": ["minimap2", "racon"],
    }.get(method, [])
    require_tools(required_exes, hint=f"--method {method} needs: {', '.join(required_exes)}. Install the missing tool(s) (e.g. `pixi add <tool>` / `conda install -c bioconda <tool>`) and make sure the correct environment is activated.")

    # the reads each method needs were checked above
    if method == "polypolish":
        assert forward_reads
        ret, out_path = run_polypolish(infile, forward_reads, reverse_reads, cpus, workdir, polish_log, debug)
    elif method == "pypolca":
        assert forward_reads
        ret, out_path = run_pypolca(infile, forward_reads, reverse_reads, cpus, memory, workdir, polish_log, polished_fasta)
    elif method == "nextpolish2":
        assert forward_reads and longreads
        ret, out_path = run_nextpolish2(infile, forward_reads, reverse_reads, longreads, cpus, workdir, polish_log, debug)
    elif method == "racon":
        assert longreads
        ret, out_path = run_racon(infile, longreads, cpus, workdir, polish_log, debug)
    else:
        raise ValueError(f"Unknown polishing method: {method}")

    # Validate the polisher's output and copy it to the requested destination
    # — done once here rather than duplicated in every run_<method>() function.
    if ret != 0 or not Path(out_path).exists() or Path(out_path).stat().st_size == 0:
        raise RuntimeError(f"{method} failed (exit {ret}); check log: {Path(workdir, polish_log)}")
    shutil.copyfile(out_path, polished_fasta)
    logger.info("AAFTF polish completed.")
    logger.info(f"{method} polished assembly: {polished_fasta}")

    next_out = next_step_name(polished_fasta, ".final.fasta")

    cleanup_workdir(workdir, debug, custom_workdir)

    if not pipe:
        logger.info(f"Your next command might be:\nAAFTF sort -i {polished_fasta} -o {next_out}")


def run_polypolish(infile: str, forward_reads: str, reverse_reads: str | None, cpus: int, workdir: str, polish_log: str, debug: bool) -> tuple[int, str]:
    """Polish with Polypolish.

    Workflow: index the assembly, align each read file separately with
    ``bwa mem -a`` (reporting all alignments, required by Polypolish), filter
    the resulting SAM pairs by insert size, then polish.

    Args:
        infile: Input assembly FASTA (copied into ``workdir``).
        forward_reads: Forward reads FASTQ.
        reverse_reads: Reverse reads FASTQ.
        cpus: Number of threads.
        workdir: Working directory.
        polish_log: Log file name (inside ``workdir``) for Polypolish stderr.
        debug: Show command output when True.

    Returns:
        Tuple of (return code, polished FASTA path). Validation, copying to the final
        destination, and status reporting are handled centrally by ``run()``.

    Raises:
        ValueError: If ``reverse_reads`` is not given.
    """
    if not reverse_reads:
        raise ValueError("--method polypolish requires paired reads (-l/--left and -r/--right)")

    assembly = str(Path(workdir, Path(infile).name))
    shutil.copyfile(infile, assembly)
    asm_name = Path(assembly).name

    bwa_index = ["bwa", "index", asm_name]
    run_cmd(bwa_index, debug, cwd=workdir)

    sam1, sam2 = "alignments_1.sam", "alignments_2.sam"
    for reads, sam_out in ((forward_reads, sam1), (reverse_reads, sam2)):
        bwa_cmd = ["bwa", "mem", "-t", str(cpus), "-a", asm_name, reads]
        with open(str(Path(workdir, sam_out)), "w") as out_fh:
            run_cmd(bwa_cmd, debug, cwd=workdir, stdout=out_fh)

    filt1, filt2 = "filtered_1.sam", "filtered_2.sam"
    filter_cmd = ["polypolish", "filter", "--in1", sam1, "--in2", sam2, "--out1", filt1, "--out2", filt2]
    run_cmd(filter_cmd, debug, cwd=workdir)

    polish_cmd = ["polypolish", "polish", asm_name, filt1, filt2]
    print_cmd(polish_cmd)
    out_path = str(Path(workdir, "polypolish_corrected.fasta"))
    with open(str(Path(workdir, polish_log)), "w") as logfile, open(out_path, "w") as out_fh:
        ret = subprocess.run(polish_cmd, cwd=workdir, stdout=out_fh, stderr=logfile)
    return ret.returncode, out_path


def run_pypolca(infile: str, forward_reads: str, reverse_reads: str | None, cpus: int, memory: int, workdir: str, polish_log: str, polished_fasta: str) -> tuple[int, str]:
    """Polish with pypolca (POLCA algorithm reimplemented in Python; runs bwa+samtools+freebayes internally).

    Copies pypolca's ``.vcf``/``.report`` sidecar outputs alongside
    ``polished_fasta`` itself (as ``{polished_fasta}.vcf`` /
    ``{polished_fasta}.pypolca_report.txt``), since those are specific to this
    method.

    Args:
        infile: Input assembly FASTA.
        forward_reads: Forward reads FASTQ.
        reverse_reads: Reverse reads FASTQ, or None.
        cpus: Number of threads.
        memory: Total memory in GB; divided by ``cpus`` for pypolca's per-thread ``-m`` (min 1G).
        workdir: Working directory.
        polish_log: Log file name (inside ``workdir``) for pypolca output.
        polished_fasta: Final output FASTA path, used to name the sidecar files.

    Returns:
        Tuple of (return code, polished FASTA path); validating and copying that primary
        output is handled centrally by ``run()``.
    """
    memperthread = f"{max(int(memory / cpus), 1)}G"  # at least 1G per thread

    pypolca_prefix = "pypolca"
    pypolca_outdir = str(Path(workdir, "pypolca_out"))
    pypolca_cmd = ["pypolca", "run", "-a", str(Path(infile).resolve()), "-1", forward_reads]
    if reverse_reads:
        pypolca_cmd.extend(["-2", reverse_reads])
    pypolca_cmd.extend(["-t", str(cpus), "-o", pypolca_outdir, "-p", pypolca_prefix, "-m", memperthread, "-f"])
    print_cmd(pypolca_cmd)
    with open(str(Path(workdir, polish_log)), "w") as logfile:
        ret = subprocess.run(pypolca_cmd, cwd=workdir, stderr=logfile, stdout=logfile)
    out_path = str(Path(pypolca_outdir, f"{pypolca_prefix}_corrected.fasta"))
    if ret.returncode == 0:
        vcf_src = str(Path(pypolca_outdir, f"{pypolca_prefix}.vcf"))
        report_src = str(Path(pypolca_outdir, f"{pypolca_prefix}.report"))
        if Path(vcf_src).exists():
            shutil.copyfile(vcf_src, f"{polished_fasta}.vcf")
        if Path(report_src).exists():
            shutil.copyfile(report_src, f"{polished_fasta}.pypolca_report.txt")
    return ret.returncode, out_path


def run_nextpolish2(infile: str, forward_reads: str, reverse_reads: str | None, longreads: str, cpus: int, workdir: str, polish_log: str, debug: bool) -> tuple[int, str]:
    """Polish with NextPolish2.

    Workflow: build a short-read k-mer (yak) database, map HiFi long reads to
    the assembly with minimap2, then run nextPolish2's repeat-aware
    correction using the HiFi alignments plus the yak k-mer database.

    Args:
        infile: Input assembly FASTA (copied into ``workdir``).
        forward_reads: Forward short reads FASTQ.
        reverse_reads: Reverse short reads FASTQ, or None.
        longreads: HiFi long-read FASTQ.
        cpus: Number of threads.
        workdir: Working directory.
        polish_log: Log file name (inside ``workdir``) for nextPolish2 output.
        debug: Show command output when True.

    Returns:
        Tuple of (return code, polished FASTA path). Validation, copying to the final
        destination, and status reporting are handled centrally by ``run()``.
    """
    assembly = str(Path(workdir, Path(infile).name))
    shutil.copyfile(infile, assembly)
    asm_name = Path(assembly).name

    yak_db = "shortreads.yak"
    yak_cmd = ["yak", "count", "-k31", "-b37", "-t", str(cpus), "-o", yak_db, forward_reads]
    if reverse_reads:
        yak_cmd.append(reverse_reads)
    run_cmd(yak_cmd, debug, cwd=workdir)

    hifi_bam = "hifi.map.bam"
    minimap_cmd = ["minimap2", "-ax", "map-hifi", "-t", str(cpus), asm_name, longreads]
    align_to_sorted_bam(minimap_cmd, str(Path(workdir, hifi_bam)), cpus, cwd=workdir, debug=debug)

    out_fasta = "nextpolish2_corrected.fasta"
    nextpolish2_cmd = ["nextPolish2", "-t", str(cpus), "-o", out_fasta, hifi_bam, asm_name, yak_db]
    print_cmd(nextpolish2_cmd)
    out_path = str(Path(workdir, out_fasta))
    with open(str(Path(workdir, polish_log)), "w") as logfile:
        ret = subprocess.run(nextpolish2_cmd, cwd=workdir, stderr=logfile, stdout=logfile)
    return ret.returncode, out_path


def run_racon(infile: str, longreads: str, cpus: int, workdir: str, polish_log: str, debug: bool) -> tuple[int, str]:
    """Polish with Racon.

    Workflow: map long reads to the assembly with minimap2 (producing PAF
    overlaps), then run racon using those overlaps to correct the assembly.

    Args:
        infile: Input assembly FASTA (copied into ``workdir``).
        longreads: Long-read FASTQ.
        cpus: Number of threads.
        workdir: Working directory.
        polish_log: Log file name (inside ``workdir``) for racon stderr.
        debug: Show command output when True.

    Returns:
        Tuple of (return code, polished FASTA path). Validation, copying to the final
        destination, and status reporting are handled centrally by ``run()``.
    """
    assembly = str(Path(workdir, Path(infile).name))
    shutil.copyfile(infile, assembly)
    asm_name = Path(assembly).name

    overlaps = "overlaps.paf"
    minimap_cmd = ["minimap2", "-x", "map-ont", "-t", str(cpus), asm_name, longreads]
    with open(str(Path(workdir, overlaps)), "w") as out_fh:
        run_cmd(minimap_cmd, debug, cwd=workdir, stdout=out_fh)

    racon_cmd = ["racon", "-t", str(cpus), longreads, overlaps, asm_name]
    print_cmd(racon_cmd)
    out_path = str(Path(workdir, "racon_corrected.fasta"))
    with open(str(Path(workdir, polish_log)), "w") as logfile, open(out_path, "w") as out_fh:
        ret = subprocess.run(racon_cmd, cwd=workdir, stdout=out_fh, stderr=logfile)
    return ret.returncode, out_path
