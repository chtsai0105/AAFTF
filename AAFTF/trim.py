"""Trims FASTQ files for reads.

This is usually for Illumina reads to quality trim reads.
This uses either fastp, includes merging step for paired reads, OR
trimmomatic. Expects adaptor sequence files to be in trimmomatic installed folder.
"""

import os
import subprocess
import sys
from os.path import dirname

from AAFTF.utility import Fzip_inplace, SafeRemove, countfastq, getRAM, printCMD, status, which_path

TRIMMOMATIC_TRUSEQSE = "adapters/TruSeq3-SE.fa"
TRIMMOMATIC_TRUSEQPE = "adapters/TruSeq3-PE.fa"
# process trimming reads with trimmomatic
# Homebrew install of trimmomatic uses a shell script
"""
#!/bin/bash -l
TRIMJAR=/usr/local/Cellar/trimmomatic/0.36/libexec/trimmomatic-0.36.jar
exec java -jar $TRIMJAR "$@"
"""
# while bioconda install uses a python script that launches java apps


def find_trimmomatic():
    """Finds the trimmomatic jar file."""
    trim_path = which_path("trimmomatic")
    if trim_path:
        with open(os.path.abspath(trim_path)) as trim_shell:
            firstLine = trim_shell.readline()
            if "#!/bin/bash" in firstLine:  # homebrew get jar location
                for line in trim_shell:
                    if line.startswith("exec java"):
                        items = line.split(" ")
                        for x in items:
                            if x.endswith(".jar"):
                                return x
            elif "#!/usr/bin/env python" in firstLine:
                trimjardir = os.path.dirname(os.path.realpath(trim_path))
                return os.path.join(trimjardir, "trimmomatic.jar")
            else:
                return False
    else:
        return False


# flake8: noqa: C901
def run(
    left,
    right=None,
    basename=None,
    method="bbduk",
    cpus=1,
    memory=None,
    minlen=75,
    avgqual=10,
    trimmomatic_adaptors="TruSeq3-PE.fa",
    trimmomatic_clip="2:30:10",
    trimmomatic_leadingwindow=3,
    trimmomatic_trailingwindow=3,
    trimmomatic_slidingwindow="4:15",
    trimmomatic_quality="phred33",
    merge=False,
    dedup=False,
    cutfront=False,
    cuttail=False,
    cutright=False,
    debug=False,
    pipe=False,
    **kwargs,
):
    """Run command for the module subtool of AAFTF."""
    if not basename:
        if "_" in os.path.basename(left):
            basename = os.path.basename(left).split("_")[0]
        elif "." in os.path.basename(left):
            basename = os.path.basename(left).split(".")[0]
        else:
            basename = os.path.basename(left)

    total = countfastq(left)
    if right:
        total = total * 2
    status(f"Loading {total:,} total reads")

    if method == "bbduk":
        if memory:
            MEM = f"-Xmx{memory}g"
        else:
            MEM = f"-Xmx{round(0.6 * getRAM())}g"

        status("Adapter trimming using BBDuk")
        bbduk_base = [
            "bbduk.sh",
            MEM,
            "ref=adapters",
            f"t={cpus}",
            "ktrim=r",
            "k=23",
            "mink=11",
            f"minlen={minlen}",
            "hdist=1",
            f"maq={avgqual}",
            "ftm=5",
            "tpe",
            "tbo",
            "overwrite=true",
        ]
        if left and right:
            # Paired mode (in1=/in2=) hits a bug in this BBDuk build's
            # PairStreamer on large/variable-length paired FASTQ: it silently
            # truncates the stream after a few hundred reads instead of
            # erroring loudly. Feed it an INTERLEAVED single file instead
            # (bbduk's single-end reader handles the full file correctly),
            # then de-interleave the trimmed output.
            interleaved_in = f"{basename}_ivl.fq.gz"
            interleaved_out = f"{basename}_ivl.trimmed.fq.gz"
            shuffle_cmd = ["shuffle.sh", f"in1={left}", f"in2={right}", f"out={interleaved_in}"]
            printCMD(shuffle_cmd)
            if debug:
                subprocess.run(shuffle_cmd)
            else:
                subprocess.run(shuffle_cmd, stderr=subprocess.DEVNULL)

            cmd = bbduk_base + [f"in={interleaved_in}", "interleaved=true", f"out={interleaved_out}"]
            printCMD(cmd)
            if debug:
                subprocess.run(cmd)
            else:
                subprocess.run(cmd, stderr=subprocess.DEVNULL)

            reformat_cmd = [
                "reformat.sh",
                f"in={interleaved_out}",
                f"out1={basename}_1P.fastq.gz",
                f"out2={basename}_2P.fastq.gz",
            ]
            printCMD(reformat_cmd)
            if debug:
                subprocess.run(reformat_cmd)
            else:
                subprocess.run(reformat_cmd, stderr=subprocess.DEVNULL)
            SafeRemove(interleaved_in)
            SafeRemove(interleaved_out)
        elif left:
            cmd = bbduk_base + [f"in={left}", f"out={basename}_1U.fastq.gz"]
            printCMD(cmd)
            if debug:
                subprocess.run(cmd)
            else:
                subprocess.run(cmd, stderr=subprocess.DEVNULL)

        if right:
            clean = countfastq(f"{basename}_1P.fastq.gz")
            clean = clean * 2
            status(f"{clean:,} reads remaining and writing to file")
            status(
                "Trimming finished:\n\tFor: {:}\n\tRev {:}".format(basename + "_1P.fastq.gz", basename + "_2P.fastq.gz")
            )
            if not pipe:
                status(
                    "Your next command might be:\n\t"
                    + "AAFTF filter -l {:} -r {:} -o {:} -c {:}\n".format(
                        basename + "_1P.fastq.gz", basename + "_2P.fastq.gz", basename, cpus
                    )
                )
        else:
            clean = countfastq(f"{basename}_1U.fastq.gz")
            status(f"{clean:,} reads remaining and writing to file")
            status("Trimming finished:\n\tSingle: {:}".format(basename + "_1U.fastq.gz"))
            if not pipe:
                status(
                    "Your next command might be:\n\t"
                    + "AAFTF filter -l {:} -o {:} -c {:}\n".format(basename + "_1U.fastq.gz", basename, cpus)
                )

    elif method == "trimmomatic":
        # find path
        trimmomatic_path = find_trimmomatic()
        if trimmomatic_path:
            jarfile = trimmomatic_path
        else:
            status("Trimmomatic cannot be found - " + "please provide location of trimmomatic.jar file.")
            sys.exit(1)

        if jarfile:
            path_to_adaptors = trimmomatic_adaptors
            clipstr = f"ILLUMINACLIP:{path_to_adaptors}:{trimmomatic_clip}"
            leadingwindow = f"LEADING:{trimmomatic_leadingwindow}"
            trailingwindow = f"TRAILING:{trimmomatic_trailingwindow}"
            slidingwindow = f"SLIDINGWINDOW:{trimmomatic_slidingwindow}"

            quality = trimmomatic_quality
            quality = f"-{quality}"  # add leading dash

            if not os.path.exists(path_to_adaptors):
                if right:
                    path_to_adaptors = os.path.join(dirname(jarfile), TRIMMOMATIC_TRUSEQPE)
                else:
                    path_to_adaptors = os.path.join(dirname(jarfile), TRIMMOMATIC_TRUSEQSE)

                if not os.path.exists(path_to_adaptors):
                    findpath = dirname(jarfile)
                    path_to_adaptors = ""
                    while findpath:
                        if os.path.exists(findpath + "/share"):
                            if right:
                                path_to_adaptors = os.path.join(findpath, "/share/trimmomatic", TRIMMOMATIC_TRUSEQPE)
                            else:
                                path_to_adaptors = os.path.join(findpath, "/share/trimmomatic", TRIMMOMATIC_TRUSEQSE)
                            break
                        new_path = dirname(findpath)
                        if new_path == findpath:  # reached filesystem root
                            break
                        findpath = new_path

                if not os.path.exists(path_to_adaptors):
                    status("Cannot find adaptors file please specify manually")
                    return

            cmd = []

            if left and right:
                cmd = [
                    "java",
                    "-jar",
                    jarfile,
                    "PE",
                    "-threads",
                    str(cpus),
                    quality,
                    left,
                    right,
                    basename + "_1P.fastq",
                    basename + "_1U.fastq",
                    basename + "_2P.fastq",
                    basename + "_2U.fastq",
                    clipstr,
                    leadingwindow,
                    trailingwindow,
                    slidingwindow,
                    f"MINLEN:{minlen}",
                ]
            elif left and not right:
                cmd = [
                    "java",
                    "-jar",
                    jarfile,
                    "SE",
                    "-threads",
                    str(cpus),
                    quality,
                    left,
                    basename + "_1U.fastq",
                    clipstr,
                    leadingwindow,
                    trailingwindow,
                    slidingwindow,
                    f"MINLEN:{minlen}",
                ]
            else:
                status("Must provide left and right pairs or single read set")
                return

            status("Running trimmomatic adapter and quality trimming")
            printCMD(cmd)
            if debug:
                subprocess.run(cmd)
            else:
                subprocess.run(cmd, stderr=subprocess.DEVNULL)
            if right:
                status("Compressing trimmed PE FASTQ files")
                Fzip_inplace(basename + "_1P.fastq", cpus)
                Fzip_inplace(basename + "_2P.fastq", cpus)
                SafeRemove(basename + "_1U.fastq")
                SafeRemove(basename + "_2U.fastq")
                status(
                    "Trimming finished:\n\tFor: {:}\n\tRev {:}".format(
                        basename + "_1P.fastq.gz", basename + "_2P.fastq.gz"
                    )
                )
                if not pipe:
                    status(
                        "Your next command might be:\n\t"
                        + "AAFTF filter -l {:} -r {:} -o {:} -c {:}\n".format(
                            basename + "_1P.fastq.gz", basename + "_2P.fastq.gz", basename, cpus
                        )
                    )
            else:
                status("Compressing trimmed SE FASTQ file")
                Fzip_inplace(basename + "_1U.fastq", cpus)
                status("Trimming finished:\n\tSingle: {:}".format(basename + "_1U.fastq.gz"))
                if not pipe:
                    status(
                        "Your next command might be:\n\t"
                        + "AAFTF filter -l {:} -o {:} -c {:}\n".format(basename + "_1U.fastq.gz", basename, cpus)
                    )

    elif method == "fastp":
        status("Adapter trimming using fastp")
        cmd = [
            "fastp",
            "--low_complexity_filter",
            "-l",
            f"{minlen}",
            "--average_qual",
            f"{avgqual}",
            "-w",
            f"{cpus}",
        ]

        #               '-wref=adapters', 't={:}'.format(cpus), 'ktrim=r',
        #           'k=23', 'mink=11', 'minlen={:}'.format(minlen), 'hdist=1',
        #           'ftm=5', 'tpe', 'tbo', 'overwrite=true']
        if left and right:
            # could add merging ...
            cmd += [
                f"--in1={left}",
                f"--in2={right}",
                f"--out1={basename}_1P.fastq.gz",
                f"--out2={basename}_2P.fastq.gz",
            ]
            if merge:
                cmd += ["--merge", f"--merged_out={basename}_MG.fastq.gz"]

        elif left:
            cmd += [f"--in={left}", f"--out={basename}_1U.fastq.gz"]
        if dedup:
            cmd += ["--dedup"]
        if cutfront:
            cmd += ["--cut_front"]
        if cuttail:
            cmd += ["--cut_tail"]
        if cutright:
            cmd += ["--cut_right"]

        cmd += [f"--html={basename}.fastp.html", f"--json={basename}.fastp.json"]
        printCMD(cmd)
        if debug:
            subprocess.run(cmd)
        else:
            subprocess.run(cmd, stderr=subprocess.DEVNULL)

        if right:
            clean = countfastq(f"{basename}_1P.fastq.gz")
            clean = clean * 2
            status(f"{clean:,} reads remaining and writing to file")
            status(
                "Trimming finished:\n\tFor: {:}\n\tRev {:}".format(basename + "_1P.fastq.gz", basename + "_2P.fastq.gz")
            )
            if not pipe:
                status(
                    "Your next command might be:\n\t"
                    + "AAFTF filter -l {:} -r {:} -o {:} -c {:}\n".format(
                        basename + "_1P.fastq.gz", basename + "_2P.fastq.gz", basename, cpus
                    )
                )
        else:
            clean = countfastq(f"{basename}_1U.fastq.gz")
            status(f"{clean:,} reads remaining and writing to file")
            status("Trimming finished:\n\tSingle: {:}".format(basename + "_1U.fastq.gz"))
            if not pipe:
                status(
                    "Your next command might be:\n\t"
                    + "AAFTF filter --left {:} -o {:} -c {:}\n".format(basename + "_1U.fastq.gz", basename, cpus)
                )

    else:
        status(f"Unknown trimming method: {method}")
