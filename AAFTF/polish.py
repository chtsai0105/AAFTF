"""Runs the pilon illumina-based polishing of assembly.

This takes care of running multiple rounds and updating the assembly
through these versions and removing temporary BAM alignment files.
"""

import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

from packaging.version import Version

from AAFTF.utility import SafeRemove, get_samtools_version, line_count, printCMD, samtools_sort_cmd, samtools_view_bam_cmd, status


def run(  # noqa: C901
    infile,
    outfile=None,
    method="pilon",
    memory=16,
    cpus=1,
    left=None,
    right=None,
    longreads=None,
    workdir=None,
    iterations=5,
    diploid=False,
    ploidy=1,
    polca="polca.sh",
    polca_samtools=None,
    debug=False,
    pipe=False,
    **kwargs,
):
    """Execute polishing step in multiple rounds provided with illumina reads and contig assembly FastA file."""
    # find reads for polishing

    polishMethod = method
    memperthread = int(memory / cpus)
    if memperthread == 0:
        # minimum should be 1Gb
        memperthread = "1"

    memperthread = f"{memperthread}G"
    status(f"calling {polishMethod} with memoryper thread as {memperthread}, input memory is {memory}GB and num cpus is {cpus}")
    forReads, revReads = (None,) * 2
    if left:
        forReads = str(Path(left).resolve())
    if right:
        revReads = str(Path(right).resolve())

    if longreads:
        # racon or nextpolish can use these
        longreads = str(Path(longreads).resolve())
    elif polishMethod == "racon":
        status(f"calling {polishMethod} without long reads (pacbio/ONT)")
        sys.exit(1)

    if method == "racon" and not longreads:
        status("Unable to located long read FASTQ raw reads, pass via -lr or --longreads")
        sys.exit(1)
    if not forReads and method != "racon":
        status("Unable to located FASTQ raw reads, pass via -l,--left and/or -r,--right")
        sys.exit(1)

    custom_workdir = 1
    if not workdir:
        custom_workdir = 0
        workdir = f"aaftf-polish_{str(uuid.uuid4())[:8]}"
    if not Path(workdir).exists():
        Path(workdir).mkdir()

    # Output file
    polishedFasta = None
    if outfile:
        polishedFasta = outfile
    else:
        fbasename = Path(infile).name.split(".f")[0]
        polishedFasta = f"{fbasename}.polished.fasta"

    method = method.lower()
    nextPolishExe = None
    polish_log = f"{method}.log"

    # Preflight: make sure the external tools this method needs are on PATH.
    # bwa is used by make_bwa_bam() for every short-read method.
    # For polca/masurca, honour a user-supplied --polca path/name rather than
    # assuming the executable is literally "polca.sh".
    polca_exe = polca
    required_exes = {
        "pilon": ["bwa", "samtools", "pilon"],
        "nextpolish": ["bwa", "samtools", "nextPolish"],
        "polca": ["bwa", "samtools", polca_exe],
        "masurca": ["bwa", "samtools", polca_exe],
    }.get(method, [])
    # shutil.which() resolves both bare names on PATH and explicit paths.
    missing = [exe for exe in required_exes if shutil.which(exe) is None]
    if missing:
        status(f"ERROR: required executable(s) not found on PATH for --method {method}: {', '.join(missing)}")
        status("Install the missing tool(s) (e.g. `pixi add pilon` / `conda install -c bioconda pilon`) and ensure the correct environment is activated.")
        sys.exit(1)

    if method == "pilon" or method == "nextpolish":
        if iterations < 1:
            status("ERROR: --iterations must be >= 1")
            sys.exit(1)
        for i in range(1, iterations + 1):
            status(f"Starting {method} polishing iteration {i}")
            correctedBase = f"polished{i}"
            if i == 1:  # first loop
                initialFasta = infile
                initialFasta = str(Path(workdir, Path(infile).name))
                shutil.copyfile(infile, initialFasta)
            else:
                initialFasta = str(Path(workdir, "polished" + str(i - 1) + ".fasta"))
            BAMfile = make_bwa_bam(initialFasta, forReads, revReads, workdir, cpus, memperthread)
            if not Path(workdir, BAMfile).exists():
                status(f"BAMfile {BAMfile} did not get created for {forReads} {revReads} in {workdir}")
                sys.exit(1)
            run_cmd = []
            dirty = []
            if method == "pilon":
                # run Pilon
                run_cmd = ["pilon", "--genome", Path(initialFasta).name, "--frags", BAMfile, f"-Xmx{memory}g", "--output", correctedBase, "--threads", str(cpus), "--changes"]
                if diploid or ploidy == 2:
                    run_cmd.append("--diploid")

                polish_log = "pilon_" + str(i) + ".log"

                printCMD(run_cmd)
                with open(str(Path(workdir, polish_log)), "w") as logfile:
                    subprocess.run(run_cmd, cwd=workdir, stderr=logfile, stdout=logfile)
                n_chg = line_count(str(Path(workdir, correctedBase + ".changes")))

                status(f"Found {n_chg:,} changes in Pilon iteration {i}")
                if n_chg == 0:
                    status(f"No changes found in Pilon iteration {i}, stopping")
                    break
            elif method == "nextpolish":
                if not nextPolishExe:
                    nextPolishmain = shutil.which("nextPolish")
                    print(nextPolishmain)
                    nextPolishExe = str(Path(Path(nextPolishmain).parent.parent, "share", "nextpolish-1.4.1", "lib", "nextpolish1.py"))
                    if not Path(nextPolishExe).exists():
                        nextPolishExe = str(Path(Path(nextPolishmain).parent, "lib", "nextpolish1.py"))
                if not nextPolishExe or not Path(nextPolishExe).exists():
                    status("Cannot find nextPolish python script")
                    return -1
                print(f"np is {nextPolishExe}")
                # initialFasta should have already been copied to the working dir
                # or is carryforward from last iteration
                run_cmd = ["samtools", "faidx", Path(initialFasta).name]
                subprocess.run(run_cmd, cwd=workdir)
                tempoutfasta = f"temp_{correctedBase}.fasta"
                run_cmd = ["python", nextPolishExe, "-g", Path(initialFasta).name, "-t", "1", "-s", BAMfile, "-p", str(cpus), "-o", tempoutfasta]
                if diploid or ploidy == 2:
                    run_cmd.extend(["-ploidy", "2"])
                elif ploidy:
                    run_cmd.extend(["-ploidy", str(ploidy)])
                polish_log = "nextpolish_t1_" + str(i) + ".log"
                with open(str(Path(workdir, polish_log)), "w") as logfile:
                    printCMD(run_cmd)
                    subprocess.run(run_cmd, cwd=workdir, stderr=logfile, stdout=logfile)
                logfile.close()
                run_cmd = ["samtools", "faidx", tempoutfasta]
                # make second BAM file for second task of nextPolish
                BAMfile = make_bwa_bam(tempoutfasta, forReads, revReads, workdir, cpus, memperthread)

                run_cmd = ["python", nextPolishExe, "-g", tempoutfasta, "-t", "2", "-debug", "-p", str(cpus), "-s", BAMfile, "-o", correctedBase + ".fasta"]
                if diploid or ploidy == 2:
                    run_cmd.extend(["-ploidy", "2"])
                elif ploidy:
                    run_cmd.extend(["-ploidy", str(ploidy)])
                polish_log = "nextpolish_t2_" + str(i) + ".log"
                with open(str(Path(workdir, polish_log)), "w") as logfile:
                    printCMD(run_cmd)
                    subprocess.run(run_cmd, cwd=workdir, stderr=logfile, stdout=logfile)
                dirty.append(tempoutfasta)

            # clean-up as we iterate to prevent tmp directory from blowing up
            dirty.extend([initialFasta + ".sa", initialFasta + ".amb", initialFasta + ".ann", initialFasta + ".pac", initialFasta + ".bwt", str(Path(workdir, BAMfile)), str(Path(workdir, BAMfile + ".bai"))])
            for f in dirty:
                if i == 1:
                    if Path(workdir, f).is_file():
                        Path(workdir, f).unlink()
                else:
                    if Path(f).is_file():
                        Path(f).unlink()

        # iteration count is the i in the counter above
        shutil.copyfile(str(Path(workdir, "polished" + str(i) + ".fasta")), polishedFasta)

        status(f"AAFTF polish completed {iterations} iterations.")
        status(f"{method} polished assembly: {polishedFasta}")
        if "_" in polishedFasta:
            nextOut = polishedFasta.split("_")[0] + ".final.fasta"
        elif "." in polishedFasta:
            nextOut = polishedFasta.split(".")[0] + ".final.fasta"
        else:
            nextOut = polishedFasta + ".final.fasta"
    elif method.lower() == "polca" or method.lower() == "masurca":
        # Warn early if the system samtools is incompatible with polca.sh
        samtoolsver = get_samtools_version()
        if samtoolsver >= Version("1.21") and not polca_samtools:
            status(f"(this may be red herring now) WARNING: samtools {samtoolsver} detected; polca.sh uses 'samtools sort -f' " "which was removed in samtools 1.21+. " "Use --polca_samtools to point to a compatible samtools (< 1.21), " "or switch to --method pilon.")

        # Build the subprocess environment, optionally injecting a compatible samtools
        polca_env = os.environ.copy()
        if polca_samtools:
            polca_samtools = str(Path(polca_samtools).resolve())
            # Most polca.sh builds honour the SAMTOOLS variable; also prepend its
            # directory to PATH as a belt-and-suspenders fallback.
            polca_env["SAMTOOLS"] = polca_samtools
            polca_env["PATH"] = str(Path(polca_samtools).parent) + ":" + polca_env.get("PATH", "")

        initialFasta = Path(infile).name
        shutil.copyfile(infile, str(Path(workdir, initialFasta)))
        polca_cmd = [polca, "-a", initialFasta, "-r", f"{forReads} {revReads}", "-t", str(cpus), "-m", memperthread]
        printCMD(polca_cmd)
        # run the polca polishing
        with open(str(Path(workdir, polish_log)), "w") as logfile:
            ret = subprocess.run(polca_cmd, cwd=workdir, stderr=logfile, stdout=logfile, env=polca_env)
        polca_out = str(Path(workdir, f"{initialFasta}.PolcaCorrected.fa"))
        if ret.returncode != 0 or not Path(polca_out).exists():
            status(f"ERROR: polca failed (exit {ret.returncode}); check log: {Path(workdir, polish_log)}")
            if samtoolsver >= Version("1.21") and not polca_samtools:
                status(f"NOTE: samtools {samtoolsver} is installed; polca.sh uses 'samtools sort -f' " "which was removed in samtools 1.21+. " "Re-run with --polca_samtools /path/to/old/samtools (< 1.21) " "or switch to --method pilon.")
            sys.exit(1)
        shutil.copyfile(polca_out, polishedFasta)
        shutil.copyfile(str(Path(workdir, f"{initialFasta}.vcf")), f"{polishedFasta}.vcf")
        shutil.copyfile(str(Path(workdir, f"{initialFasta}.report")), f"{polishedFasta}.polca_report.txt")
        status("AAFTF polish completed.")
        status(f"{method} polished assembly: {polishedFasta}")

    if "_" in polishedFasta:
        nextOut = polishedFasta.split("_")[0] + ".final.fasta"
    elif "." in polishedFasta:
        nextOut = polishedFasta.split(".")[0] + ".final.fasta"
    else:
        nextOut = polishedFasta + ".final.fasta"

    if not debug and not custom_workdir:
        SafeRemove(workdir)

    if not pipe:
        status("Your next command might be:\n" + f"\tAAFTF sort -i {polishedFasta} -o {nextOut}\n")


def make_bwa_bam(inFasta, forReads, revReads, workdir, cpus, memperthread):
    """Run BAM file generation from short reads on current assembly file to enable polishing."""
    ASMname = Path(inFasta).name
    ASMpref = Path(ASMname).stem
    BAM = ASMpref + ".bwa.bam"
    if Path(workdir, BAM).exists():
        return BAM
    tempfile_fixmate = f"{ASMpref}.fixmate.bam"
    tempfile_markdup = f"{ASMpref}.markdup.bam"
    tempfile_sort = f"{ASMpref}.sort.bam"
    tempfile_unsorted = f"{ASMpref}.unsorted.bam"
    tempfiles = [tempfile_fixmate, tempfile_markdup, tempfile_sort, tempfile_unsorted]
    bamthreads = 4
    if cpus < 4:
        bamthreads = cpus

    if not Path(workdir, BAM).is_file():
        bwa_index = ["bwa", "index", ASMname]
        printCMD(bwa_index)
        subprocess.run(bwa_index, cwd=workdir, stderr=subprocess.DEVNULL)
        bwa_cmd = ["bwa", "mem", "-t", str(cpus), ASMname, forReads]
        if revReads:
            bwa_cmd.append(revReads)

        # run BWA and pipe to samtools sort
        printCMD(bwa_cmd)
        p1 = subprocess.Popen(bwa_cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        samtoolsversion = get_samtools_version()

        if samtoolsversion < Version("1.0"):
            # run fix mate after creating BAM files from from bwa output with samtools < 1.0
            p2 = subprocess.Popen(samtools_view_bam_cmd("-", tempfile_unsorted, bamthreads), cwd=workdir, stdin=p1.stdout, stderr=subprocess.DEVNULL)
            p1.stdout.close()
            p2.communicate()
            samtools_cmd = ["samtools", "fixmate", "-r", tempfile_unsorted, tempfile_fixmate]
            subprocess.run(samtools_cmd, cwd=workdir, stderr=subprocess.DEVNULL)

            samtools_cmd = samtools_sort_cmd(tempfile_sort, tempfile_markdup, bamthreads, memory_per_thread=memperthread)
            printCMD(samtools_cmd)
            subprocess.run(samtools_cmd, cwd=workdir, stderr=subprocess.DEVNULL)
            # keep only paired reads
            samtools_cmd = samtools_view_bam_cmd(tempfile_sort, BAM, bamthreads, include_flags="0x2")
            printCMD(samtools_cmd)
            subprocess.run(samtools_cmd, cwd=workdir, stderr=subprocess.DEVNULL)

        else:
            # run fix mate directly from bwa output with samtools >= 1.0
            fixmate_fmt = "bam,level=1" if samtoolsversion >= Version("1.6") else "bam"
            samtools_cmd = ["samtools", "fixmate", "-O", fixmate_fmt, "-m", "-", tempfile_unsorted]
            p2 = subprocess.Popen(samtools_cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, stdin=p1.stdout)
            p1.stdout.close()
            p2.communicate()

            # sort to stdout, pipe into markdup
            sort_cmd = samtools_sort_cmd(tempfile_unsorted, "-", bamthreads, memory_per_thread=memperthread, tmp_prefix=ASMpref)
            printCMD(sort_cmd)
            p3 = subprocess.Popen(sort_cmd, stdout=subprocess.PIPE, cwd=workdir, stderr=subprocess.DEVNULL)
            samtools_cmd = ["samtools", "markdup", "-@", str(bamthreads), "-", tempfile_markdup]
            printCMD(samtools_cmd)
            p4 = subprocess.Popen(samtools_cmd, stdin=p3.stdout, cwd=workdir, stderr=subprocess.DEVNULL)
            p3.stdout.close()
            p4.communicate()

            # keep only paired reads
            samtools_cmd = samtools_view_bam_cmd(tempfile_markdup, BAM, bamthreads, include_flags="0x2")
            printCMD(samtools_cmd)
            subprocess.run(samtools_cmd, cwd=workdir, stderr=subprocess.DEVNULL)

        # BAM file needs to be indexed
        samtools_cmd = ["samtools", "index", "-@", str(cpus), BAM]
        printCMD(samtools_cmd)
        subprocess.run(samtools_cmd, cwd=workdir)

        for tfile in tempfiles:
            if Path(workdir, tfile).exists():
                Path(workdir, tfile).unlink()
    return BAM
