"""Utility scripts for parsing FASTA/FASTQ files and other shared helpers."""

import argparse as ap
import datetime
import gzip
import os
import re
import shutil
import subprocess
import sys
import textwrap
import urllib.request
import uuid
from itertools import islice
from pathlib import Path

import psutil
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator


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


def open_maybe_gz(path, mode="rt"):
    """Open ``path`` with gzip if its name ends in ``.gz``, else as a plain file."""
    return (gzip.open if str(path).endswith(".gz") else open)(path, mode)


def concat_files(paths, dest):
    """Concatenate files (gzipped ones are decompressed) into ``dest``."""
    with open(dest, "wb") as out:
        for path in paths:
            with open_maybe_gz(path, "rb") as fh:
                shutil.copyfileobj(fh, out)


def estimate_read_length(input):
    """Guess the read length in a FASTQ file, rounded to the nearest 10bp."""
    with open_maybe_gz(input) as infile:
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
        with open_maybe_gz(input, "rb") as fh:
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
# samtools helpers (samtools >= 1.13 required; pixi pins >= 1.24)
# ---------------------------------------------------------------------------


def samtools_sort_cmd(input_file, output_bam, threads=1, memory_per_thread=None, tmp_prefix=None, write_index=False):
    """Return a ``samtools sort`` command list.

    Args:
        input_file: Path to input SAM/BAM, or ``'-'`` to read from stdin.
        output_bam: Destination sorted BAM path.
        threads: Number of sort threads.
        memory_per_thread: Memory per thread string, e.g. ``'1G'``.
        tmp_prefix: Prefix for temporary sort files (passed as -T).
        write_index: Also write ``<output_bam>.bai`` while sorting.
    """
    cmd = ["samtools", "sort", "-@", str(threads)]
    if memory_per_thread:
        cmd += ["-m", str(memory_per_thread)]
    if write_index:
        cmd += ["--write-index", "-o", f"{output_bam}##idx##{output_bam}.bai"]
    else:
        cmd += ["-o", output_bam]
    if tmp_prefix:
        cmd += ["-T", tmp_prefix]
    cmd.append(input_file)
    return cmd


def align_to_sorted_bam(align_cmd, bam_out, threads=1, cwd=None, debug=False):
    """Pipe an aligner's SAM output straight into ``samtools sort``.

    Args:
        align_cmd: Aligner command list that writes SAM to stdout.
        bam_out: Destination sorted BAM, relative to the current directory
            (not ``cwd``).
        threads: samtools sort threads.
        cwd: Working directory for both commands.
        debug: Show both commands' stderr (hidden otherwise).

    The BAM is indexed (``<bam_out>.bai``) as it is written. Exits the program
    if either command fails, removing any partial BAM/index so a rerun does not
    mistake it for a finished alignment.
    """
    bam_path = Path(bam_out).resolve()
    sort_cmd = samtools_sort_cmd("-", str(bam_path), threads, write_index=True)
    stderr = None if debug else subprocess.DEVNULL
    printCMD(align_cmd)
    printCMD(sort_cmd)
    p1 = subprocess.Popen(align_cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=stderr)
    p2 = subprocess.Popen(sort_cmd, cwd=cwd, stdin=p1.stdout, stderr=stderr)
    p1.stdout.close()
    p2.communicate()
    if p1.wait() != 0 or p2.returncode != 0:
        bam_path.unlink(missing_ok=True)
        Path(f"{bam_path}.bai").unlink(missing_ok=True)
        status(f"ERROR: {align_cmd[0]} | samtools sort failed for {bam_out}")
        sys.exit(1)


def bam_read_count(bamfile):
    """Return (mapped, unmapped) counts of reads in a BAM file.

    Uses ``samtools flagstat`` primary counts, so secondary and supplementary
    alignments are not counted as extra reads.
    """
    stats = {}
    for line in execute(["samtools", "flagstat", "-O", "tsv", bamfile], quiet=True):
        passed, _failed, label = line.rstrip("\n").split("\t")
        stats[label] = passed
    primary, mapped = int(stats["primary"]), int(stats["primary mapped"])
    return mapped, primary - mapped


# IUPAC complement (U pairs with A); case is preserved so soft-masked bases stay lowercase.
COMPLEMENT = str.maketrans("ACGTURYKMSWBDHVNXacgturykmswbdhvnx", "TGCAAYRMKSWVBDHNXtgcaayrmkswvbdhnx")


def calc_nx(lengths, fraction=0.5):
    """Return (NX, LX) for a set of contig lengths.

    NX is the length of the contig at which the longest contigs first cover
    ``fraction`` of the total length; LX is how many contigs that takes.
    Returns (0, 0) for an empty input. The input is not modified.
    """
    target = sum(lengths) * fraction
    cumulative = 0
    for count, length in enumerate(sorted(lengths, reverse=True), 1):
        cumulative += length
        if cumulative >= target:
            return length, count
    return 0, 0


def printCMD(cmd):
    """Print out a command for debugging."""
    stringcmd = "{:}".format(" ".join(cmd))
    prefix = "\033[96mCMD:\033[00m "
    wrapper = textwrap.TextWrapper(initial_indent=prefix, width=80, subsequent_indent=" " * 8, break_long_words=False)
    print(wrapper.fill(stringcmd))


def run_cmd(cmd, debug=False, cwd=None, stdout=None, env=None, quiet_stdout=False):
    """Print a command, then run it with stderr hidden unless ``debug``.

    Args:
        cmd: Command list.
        debug: Show the command's stderr (and stdout, if ``quiet_stdout``).
        cwd: Working directory.
        stdout: stdout destination (e.g. an open file); inherited when None.
        env: Environment for the command; inherited when None.
        quiet_stdout: Also hide stdout unless ``debug``.

    Returns:
        subprocess.CompletedProcess.
    """
    printCMD(cmd)
    if quiet_stdout and not debug:
        stdout = subprocess.DEVNULL
    return subprocess.run(cmd, cwd=cwd, stdout=stdout, stderr=None if debug else subprocess.DEVNULL, env=env)


def require_tools(tools, hint=None):
    """Exit with an error naming any of ``tools`` that are not on PATH."""
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        status(f"ERROR: required tool(s) not found on PATH: {', '.join(missing)}")
        status(hint or "Install them (e.g. `conda install -c bioconda <tool>`) and make sure the correct environment is activated.")
        sys.exit(1)


def next_step_name(outfile, suffix):
    """Suggest the next step's output name: ``outfile`` up to its first '_' (else first '.') plus ``suffix``."""
    for sep in ("_", "."):
        if sep in outfile:
            return outfile.split(sep)[0] + suffix
    return outfile + suffix


def basename_from_reads(reads):
    """Derive a sample basename from a reads filename: the name up to its first '_' (else first '.')."""
    name = Path(reads).name
    for sep in ("_", "."):
        if sep in name:
            return name.split(sep)[0]
    return name


def aaftf_db_dir(required=False):
    """Return the ``$AAFTF_DB`` database directory, or None if it is unset.

    Args:
        required: Exit with an error instead of returning None when unset.
    """
    db_dir = os.environ.get("AAFTF_DB")
    if db_dir:
        return str(Path(db_dir).resolve())
    if required:
        status("ERROR: No database directory specified.\n  Set the AAFTF_DB environment variable.\n  Example:\n    export AAFTF_DB=/path/to/aaftf_db\n    AAFTF download")
        sys.exit(1)
    return None


class _Redirect308Handler(urllib.request.HTTPRedirectHandler):
    """Extend urllib's redirect handler to also follow HTTP 308.

    Python < 3.11 does not handle 308 (Permanent Redirect): the base
    ``redirect_request()`` only allows {301,302,303,307} and raises HTTPError
    for anything else, so both it and an http_error_308 dispatcher are needed.
    """

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if code == 308:
            code = 307  # method-preserving permanent redirect
        return super().redirect_request(req, fp, code, msg, headers, newurl)

    def http_error_308(self, req, fp, code, msg, headers):
        return self.http_error_302(req, fp, code, msg, headers)


URL_OPENER = urllib.request.build_opener(_Redirect308Handler())


def download_file(url, dest, force=False):
    """Download ``url`` to ``dest`` unless it already exists.

    Writes to a temporary file and renames it on success, so an interrupted
    download never leaves a partial file that a later run would reuse.

    Args:
        url: Remote URL to download.
        dest: Local file path to write.
        force: Re-download even if ``dest`` exists.

    Returns:
        ``dest``.
    """
    if Path(dest).exists() and not force:
        status(f"  Already present: {dest}")
        return dest

    status(f"  Downloading {Path(dest).name} ...")
    Path(dest).parent.mkdir(parents=True, exist_ok=True)

    tmp = f"{dest}.tmp"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "AAFTF/1.0"})
        with URL_OPENER.open(req, timeout=300) as response:
            final_url = response.geturl()
            if final_url != url:
                status(f"  Redirected to {final_url}")
            with open(tmp, "wb") as outfh:
                shutil.copyfileobj(response, outfh)
        os.replace(tmp, dest)
    except Exception as e:
        status(f"  ERROR downloading {url}: {e}")
        safe_remove(tmp)
        raise

    status(f"  Saved {dest}")
    return dest


def status(string):
    """Print out status."""
    print("\033[92m[{:}]\033[00m {:}".format(datetime.datetime.now().strftime("%b %d %I:%M %p"), string))


# from https://stackoverflow.com/questions/4417546/
# constantly-print-subprocess-output-while-process-is-running


def execute(cmd, cwd=None, debug=False, quiet=False):
    """Run a command and yield its stdout line by line.

    Args:
        cmd: Command list.
        cwd: Working directory.
        debug: Show the command's stderr (hidden otherwise).
        quiet: Don't print the command (e.g. when it runs once per contig).

    Raises:
        subprocess.CalledProcessError: If the command fails and all its output
            was read. If the caller stops reading early, the command is
            stopped and its exit status is ignored.
    """
    if not quiet:
        printCMD(cmd)
    proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, text=True, stderr=None if debug else subprocess.DEVNULL)
    finished = False
    try:
        yield from proc.stdout
        finished = True
    finally:
        proc.stdout.close()
        if not finished:
            proc.kill()
        return_code = proc.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def safe_remove(path):
    """Remove a file, symlink or directory tree; do nothing if it does not exist.

    A symlink is removed itself, never the directory it points to.
    """
    path = Path(path)
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def make_workdir(workdir, prefix):
    """Create a subcommand's working directory.

    Args:
        workdir: User-supplied ``--workdir``, or None to auto-name it
            ``aaftf-<prefix>_<random id>``.
        prefix: Subcommand name used in the auto-generated name.

    Returns:
        Tuple (workdir, custom_workdir); pass both to ``cleanup_workdir``.
    """
    custom_workdir = bool(workdir)
    if not custom_workdir:
        workdir = f"aaftf-{prefix}_{uuid.uuid4().hex[:8]}"
    Path(workdir).mkdir(parents=True, exist_ok=True)
    return workdir, custom_workdir


def cleanup_workdir(workdir, debug, custom_workdir):
    """Remove an auto-generated workdir, unless debugging or the user supplied it.

    A user-supplied ``--workdir`` (which may be shared, e.g. by ``pipeline``, or
    even the current directory) is never deleted.
    """
    if not debug and not custom_workdir:
        safe_remove(workdir)
