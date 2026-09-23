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


def fastastats(input):
    """Return (number of records, total sequence length) of a FASTA file."""
    count = length = 0
    with open(input) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
            else:
                length += len(line.rstrip())
    return count, length


def filter_fasta(fasta_in, fasta_out, keep, wrap=60):
    """Write records whose ID passes ``keep(id)`` to fasta_out, wrapped at ``wrap`` columns.

    The ID is the first whitespace-delimited word of the header (as Biopython's
    ``record.id``); the full header line is kept. The default 60-column wrap
    matches Biopython's ``SeqIO.write``.

    Returns:
        Tuple (number of records written, total length written).
    """
    count = length = 0
    with open(fasta_in) as fin, open(fasta_out, "w") as fout:
        for header, seq in SimpleFastaParser(fin):
            if keep(header.split(None, 1)[0]):
                write_fasta(fout, header, seq, wrap)
                count += 1
                length += len(seq)
    return count, length


def _count_lines(fh):
    """Count lines in a binary stream, including a final line with no newline."""
    lines = 0
    last = b"\n"
    for chunk in iter(lambda: fh.read(1 << 20), b""):
        lines += chunk.count(b"\n")
        last = chunk[-1:]
    return lines + (last != b"\n")


def countfastq(input):
    """Count the number of records in a FASTQ file (gzip or regular).

    Gzipped input is decompressed with pigz (or gzip) when available, which is
    much faster than Python's gzip module.

    Raises:
        OSError: If the file is missing or not valid gzip.
        subprocess.CalledProcessError: If pigz/gzip fails to decompress it.
    """
    decompressor = input.endswith(".gz") and (shutil.which("pigz") or shutil.which("gzip"))
    if decompressor:
        with subprocess.Popen([decompressor, "-dc", input], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as proc:
            lines = _count_lines(proc.stdout)
        if proc.returncode:
            raise subprocess.CalledProcessError(proc.returncode, [decompressor, "-dc", input])
    else:
        with (gzip.open if input.endswith(".gz") else open)(input, "rb") as fh:
            lines = _count_lines(fh)
    return lines // 4


def softwrap(string, every=60):
    """Softwrap lines in a textstring (60 columns by default, as Biopython's SeqIO.write)."""
    return "\n".join(string[i : i + every] for i in range(0, len(string), every))


def write_fasta(fh, header, seq, wrap=60):
    """Write one FASTA record to an open file handle, wrapped at ``wrap`` columns."""
    fh.write(f">{header}\n")
    if seq:
        fh.write(f"{softwrap(seq, wrap)}\n")


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
