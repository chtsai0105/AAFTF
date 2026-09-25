"""Unit tests for AAFTF/filter.py.

Covers:
  - CLI parser defaults and flag presence for the 'filter' subcommand
  - run() guard: no read1 reads → sys.exit(1)
  - bbduk command construction for paired-end and single-end reads
  - bwa/bowtie2/minimap2 command construction
  - contamdb FASTA is created from source files
  - basename auto-derivation from read1-reads filename

External network access and actual aligners are mocked throughout.
"""

import gzip
import sys
from argparse import Namespace
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from aaftf.main import main
from aaftf.resources import DATABASES, SEQ_DBS
from aaftf.utility import db_write_dir

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Helper: parse 'AAFTF filter ...' without executing the tool
# ---------------------------------------------------------------------------


def _parse_filter(argv):
    captured = {}

    def _capture(**kwargs):
        captured["args"] = Namespace(**kwargs)

    with patch.object(sys, "argv", argv):
        with patch("aaftf.filter.run", side_effect=_capture):
            main()
    return captured.get("args")


# ---------------------------------------------------------------------------
# Helper: build a minimal Namespace for filter.run()
# ---------------------------------------------------------------------------


_UNSET = object()


def _make_filter_args(tmp_path, read1=_UNSET, read2=None, aligner="bbduk", **overrides):
    if read1 is _UNSET:
        read1 = str(tmp_path / "sample_R1.fastq.gz")
    defaults = dict(
        workdir=str(tmp_path / "workdir"),
        cpus=1,
        memory=None,
        read1=read1,
        read2=read2,
        basename=None,
        aligner=aligner,
        screen_accessions=None,
        screen_urls=None,
        screen_local=None,
        debug=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


# ---------------------------------------------------------------------------
# Create stub contaminant files (valid gzip + plain) so contamdb can be built
# ---------------------------------------------------------------------------


def _write_stub_gz(path: Path, content: bytes = b">stub\nATCG\n"):
    with gzip.open(path, "wb") as f:
        f.write(content)


def _write_stub_plain(path: Path, content: bytes = b">stub\nATCG\n"):
    path.write_bytes(content)


def _mock_download(url, dest, force=False):
    """Side effect for download_file: create a stub file at *dest* and return it."""
    p = Path(dest)
    p.parent.mkdir(parents=True, exist_ok=True)
    if dest.endswith(".gz"):
        _write_stub_gz(p)
    else:
        _write_stub_plain(p)
    return dest


@pytest.fixture(autouse=True)
def _stub_databases():
    """Put stub PhiX/UniVec in the (test-isolated) database folder, as `AAFTF database` would."""
    for name in ("phix", "univec"):
        path = Path(db_write_dir(), DATABASES[name]["filename"])
        _write_stub_gz(path) if path.name.endswith(".gz") else _write_stub_plain(path)


# ---------------------------------------------------------------------------
# Parser defaults and required flags
# ---------------------------------------------------------------------------


class TestFilterParser:
    def test_debug_false_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.debug is False

    def test_debug_flag_sets_true(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq", "-v"])
        assert args.debug is True

    def test_default_aligner_is_bbduk(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.aligner == "bbduk"

    def test_aligner_bowtie2(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq", "--aligner", "bowtie2"])
        assert args.aligner == "bowtie2"

    def test_aligner_bwa(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq", "--aligner", "bwa"])
        assert args.aligner == "bwa"

    def test_aligner_minimap2(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq", "--aligner", "minimap2"])
        assert args.aligner == "minimap2"

    def test_default_cpus(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.cpus == 1

    def test_custom_cpus(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq", "-c", "4"])
        assert args.cpus == 4

    def test_screen_accessions_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.screen_accessions is None

    def test_screen_urls_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.screen_urls is None

    def test_screen_local_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.screen_local is None

    def test_parses_read1_reads(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.read1 == "R1.fq"

    def test_parses_read2_reads(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq", "-2", "R2.fq"])
        assert args.read2 == "R2.fq"

    def test_read2_none_by_default(self):
        args = _parse_filter(["AAFTF", "filter", "-1", "R1.fq"])
        assert args.read2 is None

    def test_missing_read1_exits_nonzero(self):
        with patch.object(sys, "argv", ["AAFTF", "filter"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code != 0

    def test_filter_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "filter", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0


# ---------------------------------------------------------------------------
# run() guard: no read1 reads
# ---------------------------------------------------------------------------


class TestFilterRunGuards:
    def test_no_read1_raises(self, tmp_path):
        args = _make_filter_args(tmp_path, read1=None)
        from aaftf.filter import run

        with patch("aaftf.filter.download_file", side_effect=_mock_download):
            with patch("aaftf.filter.count_fastq", return_value=0):
                with pytest.raises(ValueError):
                    run(**vars(args))


# ---------------------------------------------------------------------------
# contamdb creation
# ---------------------------------------------------------------------------


class TestFilterContamdbCreation:
    def test_contamdb_created_in_workdir(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        args = _make_filter_args(tmp_path, read1=read1, aligner="bbduk")
        workdir = Path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)

        from aaftf.filter import run

        with patch("aaftf.filter.download_file", side_effect=_mock_download):
            with patch("aaftf.filter.count_fastq", return_value=100):
                with patch("aaftf.utility.subprocess.run"):
                    run(**vars(args))

        contamdb = workdir / "contamdb.fa"
        assert contamdb.exists()

    def test_screen_local_added_to_contamdb(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        local_fa = tmp_path / "extra.fa"
        local_fa.write_text(">extra\nATCGATCG\n")
        args = _make_filter_args(tmp_path, read1=read1, aligner="bbduk", screen_local=[str(local_fa)])
        Path(args.workdir).mkdir(parents=True, exist_ok=True)

        from aaftf.filter import run

        with patch("aaftf.filter.download_file", side_effect=_mock_download):
            with patch("aaftf.filter.count_fastq", return_value=100):
                with patch("aaftf.utility.subprocess.run"):
                    run(**vars(args))

        contamdb = Path(args.workdir) / "contamdb.fa"
        content = contamdb.read_text()
        assert ">extra" in content

    def test_contamdb_contains_downloaded_stubs(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        args = _make_filter_args(tmp_path, read1=read1, aligner="bbduk")
        Path(args.workdir).mkdir(parents=True, exist_ok=True)

        from aaftf.filter import run

        with patch("aaftf.filter.download_file", side_effect=_mock_download):
            with patch("aaftf.filter.count_fastq", return_value=100):
                with patch("aaftf.utility.subprocess.run"):
                    run(**vars(args))

        contamdb = Path(args.workdir) / "contamdb.fa"
        assert contamdb.stat().st_size > 0


# ---------------------------------------------------------------------------
# bbduk command construction
# ---------------------------------------------------------------------------


def _run_filter_bbduk(tmp_path, read1, read2=None, **extra):
    """Run filter.run() with aligner=bbduk; return captured subprocess commands."""
    args = _make_filter_args(tmp_path, read1=read1, read2=read2, aligner="bbduk", **extra)
    Path(args.workdir).mkdir(parents=True, exist_ok=True)
    cmds = []

    from aaftf.filter import run

    with patch("aaftf.filter.download_file", side_effect=_mock_download):
        with patch("aaftf.filter.count_fastq", return_value=100):
            with patch("aaftf.utility.subprocess.run", side_effect=lambda cmd, **kw: cmds.append(cmd)):
                run(**vars(args))
    return cmds, args


class TestFilterRunBbduk:
    # PE bbduk runs as shuffle.sh -> bbduk.sh -> reformat.sh (see a136d93,
    # the BBDuk PairStreamer paired-mode workaround), so a given argument may
    # land on any of the three commands rather than the first.
    def test_pe_command_includes_in(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1, read2)
        assert any(f"in1={read1}" in " ".join(c) for c in cmds)

    def test_pe_command_includes_in2(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1, read2)
        assert any(f"in2={read2}" in " ".join(c) for c in cmds)

    def test_pe_output_files_use_filtered_basename(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1, read2)
        assert any("sample_filtered_1.fastq.gz" in " ".join(c) for c in cmds)

    def test_se_output_uses_u_suffix(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1, read2=None)
        cmd_str = " ".join(cmds[0])
        assert "sample_filtered_U.fastq.gz" in cmd_str

    def test_command_starts_with_bbduk(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1, read2)
        assert any(c[0] == "bbduk.sh" for c in cmds)

    def test_command_includes_contamdb_ref(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, args = _run_filter_bbduk(tmp_path, read1, read2)
        # contamdb.fa is passed via ref= argument on the bbduk.sh step
        assert any("contamdb.fa" in " ".join(c) for c in cmds)

    def test_basename_derived_from_underscore_split(self, tmp_path):
        read1 = str(tmp_path / "MySample_R1.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1)
        assert any("MySample_filtered_U.fastq.gz" in " ".join(c) for c in cmds)

    def test_explicit_basename_preserved(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        cmds, _ = _run_filter_bbduk(tmp_path, read1, basename="custom")
        assert any("custom_filtered_U.fastq.gz" in " ".join(c) for c in cmds)


# ---------------------------------------------------------------------------
# bwa command construction
# ---------------------------------------------------------------------------


def _run_filter_bwa(tmp_path, read1, read2=None, **extra):
    """Run filter.run() with aligner=bwa; return captured subprocess commands."""
    args = _make_filter_args(tmp_path, read1=read1, read2=read2, aligner="bwa", **extra)
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    cmds = []
    popen_cmds = []

    mock_proc = MagicMock()
    mock_proc.stdout = MagicMock()
    mock_proc.communicate.return_value = (b"", b"")
    mock_proc.wait.return_value = 0
    mock_proc.returncode = 0

    def _fake_run(cmd, **kw):
        cmds.append(cmd)
        # When samtools sort is called, create the alignBAM so post-processing triggers
        if cmd and len(cmd) > 1 and "samtools" in str(cmd[0]):
            # find the output BAM argument (samtools sort -o <out>)
            for i, tok in enumerate(cmd):
                if tok == "-o" and i + 1 < len(cmd):
                    Path(cmd[i + 1]).touch()
                    break
            # fallback: touch any .bam in workdir
            for f in workdir.glob("*.bam"):
                pass  # already exists via above
            else:
                # create a stub alignBAM so isfile check passes
                (workdir / "_stub.bam").touch()

    from aaftf.filter import run

    with patch("aaftf.filter.download_file", side_effect=_mock_download):
        with patch("aaftf.filter.count_fastq", return_value=100):
            with patch("aaftf.utility.subprocess.run", side_effect=_fake_run):
                with patch("aaftf.utility.subprocess.Popen", side_effect=lambda cmd, **kw: (popen_cmds.append(cmd), mock_proc)[1]):
                    with patch("aaftf.filter.bam_read_count", return_value=(50, 50)):
                        with patch("aaftf.utility.safe_remove"):
                            run(**vars(args))
    return cmds, popen_cmds, args


class TestFilterRunBwa:
    def test_bwa_index_called(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, popen_cmds, _ = _run_filter_bwa(tmp_path, read1, read2)
        all_cmds = cmds + popen_cmds
        assert any(c and c[0] == "bwa" and "index" in c for c in all_cmds)

    def test_bwa_mem_called(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, popen_cmds, _ = _run_filter_bwa(tmp_path, read1, read2)
        all_cmds = cmds + popen_cmds
        assert any(c and c[0] == "bwa" and "mem" in c for c in all_cmds)

    def test_bwa_mem_includes_reads(self, tmp_path):
        read1 = str(tmp_path / "sample_R1.fastq.gz")
        read2 = str(tmp_path / "sample_R2.fastq.gz")
        cmds, popen_cmds, _ = _run_filter_bwa(tmp_path, read1, read2)
        mem_cmds = [c for c in popen_cmds if c and c[0] == "bwa" and "mem" in c]
        assert len(mem_cmds) > 0
        assert read1 in mem_cmds[0]


class TestFilterDatabases:
    def _run(self, tmp_path, **extra):
        args = _make_filter_args(tmp_path, read1=str(tmp_path / "sample_R1.fastq.gz"), **extra)
        from aaftf.filter import run

        calls = []
        with patch("aaftf.filter.download_file", side_effect=lambda url, dest, force=False: (calls.append(url), _mock_download(url, dest))[1]):
            with patch("aaftf.filter.count_fastq", return_value=100):
                with patch("aaftf.utility.subprocess.run"):
                    run(**vars(args))
        return calls

    def test_missing_phix_raises_with_download_command(self, tmp_path):
        Path(db_write_dir(), DATABASES["phix"]["filename"]).unlink()
        with pytest.raises(FileNotFoundError) as exc:
            self._run(tmp_path)
        assert "AAFTF database phix" in str(exc.value)

    def test_databases_are_not_downloaded_by_filter(self, tmp_path):
        assert self._run(tmp_path) == []

    def test_screen_accession_fetched_from_eutils(self, tmp_path):
        calls = self._run(tmp_path, screen_accessions=["NC_001422"])
        assert calls == [SEQ_DBS["nucleotide"] % "NC_001422"]
        assert calls[0].startswith("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?")
