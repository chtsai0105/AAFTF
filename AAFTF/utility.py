"""Utility scripts for parsing FASTA/FASTQ files and other shared helpers."""

import argparse as ap
import datetime
import gzip
import os
import re
import shutil
import subprocess
import textwrap
from itertools import islice
from pathlib import Path

import psutil
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from packaging.version import Version


class CustomHelpFormatter(ap.HelpFormatter):
    """Custom help formatter for argparse to enhance text wrapping and default value display.

    This formatter adjusts the text wrapping for better readability and ensures that the default value of an argument is displayed
    in the help message if not already present.
    """

    def _fill_text(self, text, width, indent):
        """Format the class/function docstring with appropriate wrapping."""
        text = [self._whitespace_matcher.sub(" ", paragraph.strip()) for paragraph in text.split("\n\n") if paragraph.strip()]
        return "\n\n".join([textwrap.fill(line, width) for line in text])

    def _split_lines(self, text, width):
        """Enable multi-line display in argument help message."""
        text = [self._whitespace_matcher.sub(" ", line.strip()) for line in text.split("\n") if line.strip()]
        return [wrapped_line for line in text for wrapped_line in textwrap.wrap(line, width)]

    def _get_help_string(self, action):
        """Allow additional message after default parameter displayed."""
        help = action.help
        pattern = r"\(default: .+\)"
        if re.search(pattern, action.help) is None:
            if action.default is not ap.SUPPRESS and action.default is not None and action.default is not False:
                defaulting_nargs = [ap.OPTIONAL, ap.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += " (default: %(default)s)"
        return help


def estimate_read_length(input):
    """Guess the read length in a FASTQ file, rounded to the nearest 10bp."""
    opener = gzip.open if input.endswith(".gz") else open
    with opener(input, "rt") as infile:
        records = islice(FastqGeneralIterator(infile), 500)
        max_len = max(len(seq) for _, seq, _ in records)
    return round(max_len, -1)


def checkfile(input):
    """Check that file to read is valid."""

    def _getSize(filename):
        st = os.stat(filename)
        return st.st_size

    if Path(input).is_file():
        return _getSize(input) >= 1
    return False


def getRAM(max_lim=0):
    """Get the available RAM on system, in GB, kept safely under the true value.

    Rounded down to the nearest 10 once available RAM exceeds 10 GB;
    otherwise rounded down to 1 decimal with a small safety margin
    subtracted, clamped at 0. This never exceeds the true available RAM.

    This is a simplistic approach which does not take into account shared
    HPC allocations.

    Args:
        max_lim: Optional upper bound (GB) to cap the result at.
    """
    avail_gb = psutil.virtual_memory().available / (1024.0**3)
    safe_gb = avail_gb // 10 * 10 if avail_gb > 10 else max(round(avail_gb, 1) - 0.1, 0)
    return min(safe_gb, max_lim) if max_lim else safe_gb


def countfasta(input):
    """Count the number of records in FastA file."""
    count = 0
    with open(input) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def fastastats(input):
    """Calculate statistics (num and total length) of FastA file."""
    count = 0
    length = 0
    with open(input) as f:
        for Header, Seq in SimpleFastaParser(f):
            count += 1
            length += len(Seq)
    return count, length


def countfastq(input):
    """Count the number of records in a FASTQ file (gzip or regular)."""
    opener = gzip.open if input.endswith(".gz") else open
    with opener(input, "rt") as fh:
        lines = sum(1 for _ in fh)
    return lines // 4


def softwrap(string, every=80):
    """Softwrap lines in a textstring."""
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i : i + every])
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# samtools version detection and version-aware command builders
# ---------------------------------------------------------------------------

_SAMTOOLS_VERSION_CACHE = None


def get_samtools_version():
    """Return the installed samtools version as a packaging.version.Version.

    The result is cached after the first call so that subsequent calls within
    the same process incur no subprocess overhead.  Returns Version("0.0") if
    samtools is not found or the version string cannot be parsed.
    """
    global _SAMTOOLS_VERSION_CACHE
    if _SAMTOOLS_VERSION_CACHE is not None:
        return _SAMTOOLS_VERSION_CACHE
    try:
        result = subprocess.run(["samtools"], capture_output=True, text=True)
        m = re.search(r"Version:\s+(\S+)", result.stderr or "")
    except (FileNotFoundError, OSError):
        m = None
    _SAMTOOLS_VERSION_CACHE = Version(m.group(1)) if m else Version("0.0")
    return _SAMTOOLS_VERSION_CACHE


def samtools_sort_cmd(input_file, output_bam, threads=1, memory_per_thread=None, tmp_prefix=None):
    """Return a samtools sort command list compatible with the installed version.

    Key version boundaries:
      < 1.0  : no -@ threads flag; sort uses positional ``infile outprefix``
      1.0-1.2: -@ threads added; output still positional prefix
      >= 1.3 : -o outfile flag available (replaces positional prefix)
      >= 1.21: -f flag removed — code that used -f breaks on this version

    Args:
        input_file: Path to input SAM/BAM, or ``'-'`` to read from stdin.
        output_bam: Destination sorted BAM path.
        threads: Number of sort threads.
        memory_per_thread: Memory per thread string, e.g. ``'1G'``.
        tmp_prefix: Prefix for temporary sort files (passed as -T).

    Returns:
        List of strings suitable for subprocess.run / subprocess.Popen.
    """
    ver = get_samtools_version()
    cmd = ["samtools", "sort"]
    if ver >= Version("1.0"):
        cmd += ["-@", str(threads)]
    if memory_per_thread:
        cmd += ["-m", str(memory_per_thread)]
    if ver >= Version("1.3"):
        # Modern syntax: samtools sort [-@ N] [-m mem] [-T pfx] -o out.bam in
        cmd += ["-o", output_bam]
        if tmp_prefix:
            cmd += ["-T", tmp_prefix]
        cmd.append(input_file)
    else:
        # Pre-1.3 (and the -f flag removed in 1.21 lives only in this branch,
        # which is unreachable on 1.21): samtools sort [opts] -f in.bam out.bam
        cmd += ["-f", input_file, output_bam]
    return cmd


def samtools_view_bam_cmd(input_file, output_bam, threads=1, include_flags=None):
    """Return a samtools view command list that produces a sorted BAM file.

    Key version boundaries:
      < 1.0  : no -@ flag; needs explicit -S to declare SAM input
      1.0-1.5: -@ flag; -b for BAM output; -S ignored (auto-detect format)
      >= 1.6 : -O bam,level=1 for compressed BAM output

    Args:
        input_file: Path to input file, or ``'-'`` for stdin.
        output_bam: Destination BAM path.
        threads: Number of threads.
        include_flags: Optional SAM flag mask for ``-f`` (integer or string).

    Returns:
        List of strings suitable for subprocess.run / subprocess.Popen.
    """
    ver = get_samtools_version()
    if ver >= Version("1.6"):
        cmd = ["samtools", "view", "-O", "bam,level=1", "-@", str(threads), "-o", output_bam]
    elif ver >= Version("1.0"):
        cmd = ["samtools", "view", "-b", "-@", str(threads), "-o", output_bam]
    else:
        # Pre-1.0: no threads flag; need -S to declare SAM input
        cmd = ["samtools", "view", "-bS", "-o", output_bam]
    if include_flags is not None:
        cmd += ["-f", str(include_flags)]
    cmd.append(input_file)
    return cmd


def bam_read_count(bamfile):
    """Count the number of reads in a BAM file using samtools."""
    cmd = ["samtools", "idxstats", bamfile]
    mapped = 0
    unmapped = 0
    for line in execute(cmd, "."):
        rname, rlen, nm, nu = line.rstrip().split()
        mapped += int(nm)
        unmapped += int(nu)
    return (mapped, unmapped)


def RevComp(s):
    """Reverse complement a DNA string."""
    rev_comp_lib = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A", "M": "K", "R": "Y", "W": "W", "S": "S", "Y": "R", "K": "M", "V": "B", "H": "D", "D": "H", "B": "V", "X": "X", "N": "N"}
    cseq = ""
    n = len(s)
    s = s.upper()
    for i in range(0, n):
        c = s[n - i - 1]
        cseq += rev_comp_lib[c]
    return cseq


def calcN50(lengths, num=0.5):
    """Calculate the N50 from a set of integers."""
    lengths.sort()
    total_len = sum(lengths)
    n50 = 0
    cumulsum = 0
    for n in reversed(lengths):
        cumulsum += n
        if n50 == 0 and cumulsum >= total_len * num:
            n50 = n
    return n50


def printCMD(cmd):
    """Print out a command for debugging."""
    stringcmd = "{:}".format(" ".join(cmd))
    prefix = "\033[96mCMD:\033[00m "
    wrapper = textwrap.TextWrapper(initial_indent=prefix, width=80, subsequent_indent=" " * 8, break_long_words=False)
    print(wrapper.fill(stringcmd))


def status(string):
    """Print out status."""
    print("\033[92m[{:}]\033[00m {:}".format(datetime.datetime.now().strftime("%b %d %I:%M %p"), string))


# from https://stackoverflow.com/questions/4417546/
# constantly-print-subprocess-output-while-process-is-running


def execute(cmd, dir):
    """Execute a command and wait for result."""
    popen = subprocess.Popen(cmd, cwd=dir, stdout=subprocess.PIPE, universal_newlines=True, stderr=subprocess.DEVNULL)
    yield from iter(popen.stdout.readline, "")
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def Fzip_inplace(input, cpus):
    """Function to run zip as fast as it can, pigz -> gzip."""
    if shutil.which("pigz"):
        cmd = ["pigz", "-f", "-p", str(cpus), input]
    else:
        cmd = ["gzip", "-f", input]
    try:
        runSubprocess(cmd, ".", log)
    except NameError:
        subprocess.call(cmd)


def SafeRemove(input):
    """Test and remove a folder or file."""
    if Path(input).is_dir():
        shutil.rmtree(input)
    elif Path(input).is_file():
        Path(input).unlink()
    else:
        return


def open_pipe(command, mode="r", buff=1024 * 1024):
    """Open read or write pipe to program.

    `command` may be an argv list (run directly, no shell) or a string
    (run via the shell) -- the string form exists for `zopen()`'s "!"
    passthrough where the caller explicitly supplies a shell command.
    """
    import signal
    import subprocess

    shell = isinstance(command, str)
    if "r" in mode:
        return subprocess.Popen(command, shell=shell, bufsize=buff, stdout=subprocess.PIPE, preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL)).stdout
    elif "w" in mode:
        return subprocess.Popen(command, shell=shell, bufsize=buff, stdin=subprocess.PIPE).stdin
    return None


NORMAL = 0
PROCESS = 1
PARALLEL = 2
WHICH_GZIP = shutil.which("gzip")
WHICH_PIGZ = shutil.which("pigz")

# streaming parallel pigz open via
# https://github.com/DarkoVeberic/utl/blob/master/futile/futile.py


def open_gz(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    """Open a gzip file using processes (gzip/pigz) or native library."""
    import gzip

    if "w" in mode:
        return gzip.open(filename, mode)

    if external is None or external == NORMAL:
        return gzip.GzipFile(filename, mode, buff)
    elif external == PROCESS:
        if not WHICH_GZIP:
            return open_gz(filename, mode, buff, NORMAL)
        if "r" in mode:
            return open_pipe([WHICH_GZIP, "-dc", filename], mode, buff)
    elif external == PARALLEL:
        if not WHICH_PIGZ:
            return open_gz(filename, mode, buff, PROCESS)
        if "r" in mode:
            return open_pipe([WHICH_PIGZ, "-dc", filename], mode, buff)
    return None


def zopen(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    """Open pipe, zipped, or unzipped file automagically.

    external == 0: normal zip libraries
    external == 1: (zcat, gzip) or (bzcat, bzip2)
    external == 2: (pigz -dc, pigz) or (pbzip2 -dc, pbzip2)
    """
    if "r" in mode and "w" in mode:
        return None
    if filename.startswith("!"):
        return open_pipe(filename[1:], mode, buff)
    elif filename.endswith(".gz"):
        return open_gz(filename, mode, buff, external)
    else:
        return open(filename, mode, buff)
    return None
