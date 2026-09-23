"""Unit tests for AAFTF/utility.py.

All tests are pure Python — no external bioinformatics tools required.
"""

import gzip
import io
import shutil
import subprocess

import pytest

from AAFTF.utility import (
    RevComp,
    SafeRemove,
    align_to_sorted_bam,
    bam_read_count,
    calcN50,
    checkfile,
    countfastq,
    fastastats,
    filter_fasta,
    samtools_sort_cmd,
    softwrap,
    write_fasta,
)
from tests.conftest import make_fastq_text

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# calcN50
# ---------------------------------------------------------------------------


class TestCalcN50:
    """calcN50(lengths) uses the default num=0.5 for N50."""

    # lengths [20, 80, 150] → total=250, 50 %=125
    # cumsum from largest: 150 ≥ 125 → N50=150
    def test_n50_basic(self):
        assert calcN50([20, 80, 150]) == 150

    # 90 %=225 → cumsum 150 < 225, 150+80=230 ≥ 225 → N90=80
    def test_n90(self):
        assert calcN50([20, 80, 150], num=0.9) == 80

    def test_single_contig(self):
        assert calcN50([500]) == 500

    def test_equal_contigs(self):
        # [100]*4 → total=400, 50 %=200; cumsum: 100 < 200, 200 ≥ 200 → N50=100
        assert calcN50([100, 100, 100, 100]) == 100

    def test_two_contigs(self):
        # [300, 100] → total=400, 50%=200; cumsum: 300 ≥ 200 → N50=300
        assert calcN50([300, 100]) == 300

    def test_input_order_irrelevant(self):
        # Function sorts internally
        assert calcN50([150, 20, 80]) == calcN50([20, 80, 150])


# ---------------------------------------------------------------------------
# RevComp
# ---------------------------------------------------------------------------


class TestRevComp:
    def test_simple(self):
        assert RevComp("ATCG") == "CGAT"

    def test_complement_only(self):
        assert RevComp("AAAA") == "TTTT"
        assert RevComp("CCCC") == "GGGG"

    def test_palindrome(self):
        assert RevComp("AATTAATT") == "AATTAATT"

    def test_uppercase_conversion(self):
        # RevComp uppercases input before processing
        assert RevComp("atcg") == "CGAT"

    def test_longer_sequence(self):
        assert RevComp("ATCGATCG") == "CGATCGAT"

    def test_all_bases(self):
        # A↔T, C↔G
        assert RevComp("ACGT") == "ACGT"


# ---------------------------------------------------------------------------
# softwrap
# ---------------------------------------------------------------------------


class TestSoftwrap:
    def test_short_string_unchanged(self):
        assert softwrap("ATCG", 80) == "ATCG"

    def test_wraps_at_specified_width(self):
        result = softwrap("A" * 20, 10)
        lines = result.split("\n")
        assert len(lines) == 2
        assert all(len(line) == 10 for line in lines)

    def test_exact_width_no_wrap(self):
        assert softwrap("A" * 10, 10) == "A" * 10

    def test_one_over_width_wraps(self):
        result = softwrap("A" * 11, 10)
        assert "\n" in result
        assert result.split("\n")[0] == "A" * 10

    def test_empty_string(self):
        assert softwrap("", 80) == ""

    def test_default_width_is_60(self):
        assert softwrap("A" * 61).split("\n") == ["A" * 60, "A"]


class TestWriteFasta:
    def test_wraps_at_60(self):
        buf = io.StringIO()
        write_fasta(buf, "ctg1 desc", "A" * 61)
        assert buf.getvalue() == ">ctg1 desc\n" + "A" * 60 + "\nA\n"

    def test_empty_sequence_has_no_blank_line(self):
        buf = io.StringIO()
        write_fasta(buf, "empty", "")
        assert buf.getvalue() == ">empty\n"


# ---------------------------------------------------------------------------
# checkfile
# ---------------------------------------------------------------------------


class TestCheckfile:
    def test_existing_nonempty_file_is_true(self, fasta_file):
        assert checkfile(str(fasta_file)) is True

    def test_empty_file_is_false(self, empty_file):
        assert checkfile(str(empty_file)) is False

    def test_missing_file_is_false(self, tmp_path):
        assert checkfile(str(tmp_path / "nonexistent.fa")) is False

    def test_symlink_to_nonempty_file_is_true(self, tmp_path, fasta_file):
        link = tmp_path / "link.fa"
        link.symlink_to(fasta_file)
        assert checkfile(str(link)) is True

    def test_broken_symlink_is_false(self, tmp_path):
        link = tmp_path / "broken_link.fa"
        link.symlink_to(tmp_path / "nonexistent_target.fa")
        assert checkfile(str(link)) is False


# ---------------------------------------------------------------------------
# fastastats / filter_fasta
# ---------------------------------------------------------------------------


class TestFastastats:
    def test_count(self, fasta_file):
        count, total_len = fastastats(str(fasta_file))
        assert count == 3

    def test_total_length(self, fasta_file):
        # SEQ1=20, SEQ2=80, SEQ3=150 → 250
        count, total_len = fastastats(str(fasta_file))
        assert total_len == 250

    def test_wrapped_sequence_length(self, tmp_path):
        p = tmp_path / "wrapped.fa"
        p.write_text(">a desc\nACGT\nAC\n>b\nGG\n")
        assert fastastats(str(p)) == (2, 8)


class TestFilterFasta:
    def test_drops_by_id_and_returns_stats(self, tmp_path):
        src = tmp_path / "in.fa"
        src.write_text(">a desc\nACGT\nAC\n>b\nGG\n>c\nTTT\n")
        out = tmp_path / "out.fa"
        assert filter_fasta(str(src), str(out), lambda seq_id: seq_id != "b") == (2, 9)
        assert out.read_text() == ">a desc\nACGTAC\n>c\nTTT\n"

    def test_matches_biopython_seqio_write(self, tmp_path):
        from Bio import SeqIO

        src = tmp_path / "in.fa"
        src.write_text(">a desc\n" + "ACGT" * 40 + "\n>b\n" + "G" * 60 + "\n>c\n" + "T" * 61 + "\n")
        out = tmp_path / "out.fa"
        expected = tmp_path / "expected.fa"
        filter_fasta(str(src), str(out), lambda seq_id: True)
        SeqIO.write(SeqIO.parse(str(src), "fasta"), str(expected), "fasta")
        assert out.read_text() == expected.read_text()

    def test_keep_none(self, tmp_path):
        src = tmp_path / "in.fa"
        src.write_text(">a\nACGT\n")
        out = tmp_path / "out.fa"
        assert filter_fasta(str(src), str(out), lambda seq_id: False) == (0, 0)
        assert out.read_text() == ""


# ---------------------------------------------------------------------------
# countfastq
# ---------------------------------------------------------------------------


class TestCountFastq:
    def test_plain_fastq(self, fastq_file):
        assert countfastq(str(fastq_file)) == 10

    def test_gzip_fastq(self, gz_fastq_file):
        assert countfastq(str(gz_fastq_file)) == 10

    def test_single_read(self, tmp_path):
        p = tmp_path / "one.fastq"
        p.write_text(make_fastq_text(1))
        assert countfastq(str(p)) == 1

    def test_gz_50_reads(self, tmp_path):
        p = tmp_path / "fifty.fastq.gz"
        with gzip.open(p, "wt") as fh:
            fh.write(make_fastq_text(50))
        assert countfastq(str(p)) == 50

    def test_missing_trailing_newline(self, tmp_path):
        p = tmp_path / "nonl.fastq"
        p.write_text(make_fastq_text(3).rstrip("\n"))
        assert countfastq(str(p)) == 3

    def test_gz_without_pigz_or_gzip(self, gz_fastq_file, monkeypatch):
        monkeypatch.setattr("AAFTF.utility.shutil.which", lambda name: None)
        assert countfastq(str(gz_fastq_file)) == 10

    def test_bad_gzip_raises(self, tmp_path):
        p = tmp_path / "bad.fastq.gz"
        p.write_bytes(b"not a valid gzip file")
        with pytest.raises((OSError, subprocess.CalledProcessError)):
            countfastq(str(p))


# ---------------------------------------------------------------------------
# SafeRemove
# ---------------------------------------------------------------------------


class TestSafeRemove:
    def test_remove_file(self, tmp_path):
        p = tmp_path / "todelete.txt"
        p.write_text("bye")
        SafeRemove(str(p))
        assert not p.exists()

    def test_remove_directory(self, tmp_path):
        d = tmp_path / "subdir"
        d.mkdir()
        (d / "file.txt").write_text("hi")
        SafeRemove(str(d))
        assert not d.exists()

    def test_remove_nonexistent_is_noop(self, tmp_path):
        # Should not raise
        SafeRemove(str(tmp_path / "ghost"))


# ---------------------------------------------------------------------------
# samtools helpers
# ---------------------------------------------------------------------------


class TestSamtoolsSortCmd:
    def test_basic(self):
        assert samtools_sort_cmd("-", "out.bam", 4) == ["samtools", "sort", "-@", "4", "-o", "out.bam", "-"]

    def test_memory_and_tmp_prefix(self):
        cmd = samtools_sort_cmd("in.sam", "out.bam", 2, memory_per_thread="1G", tmp_prefix="tmp")
        assert cmd == ["samtools", "sort", "-@", "2", "-m", "1G", "-o", "out.bam", "-T", "tmp", "in.sam"]

    def test_write_index_requests_bai(self):
        cmd = samtools_sort_cmd("-", "out.bam", 1, write_index=True)
        assert cmd == ["samtools", "sort", "-@", "1", "--write-index", "-o", "out.bam##idx##out.bam.bai", "-"]


_SAM = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:100\nr1\t0\tc1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
# 2 reads: r1 mapped (plus a supplementary and a secondary record), r2 unmapped
_SAM_MIXED = _SAM + "r1\t2048\tc1\t50\t60\t4M\t*\t0\t0\tACGT\tIIII\nr1\t256\tc1\t30\t0\t4M\t*\t0\t0\tACGT\tIIII\nr2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not installed")
class TestAlignToSortedBam:
    def test_writes_bam_relative_to_current_dir_not_cwd(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        sub = tmp_path / "work"
        sub.mkdir()
        (sub / "in.sam").write_text(_SAM)
        align_to_sorted_bam(["cat", "in.sam"], "out.bam", cwd=str(sub), stderr=subprocess.DEVNULL)
        assert (tmp_path / "out.bam").stat().st_size > 0
        assert (tmp_path / "out.bam.bai").exists()
        assert not (sub / "out.bam").exists()

    def test_bam_read_count_ignores_secondary_and_supplementary(self, tmp_path):
        (tmp_path / "in.sam").write_text(_SAM_MIXED)
        bam = tmp_path / "out.bam"
        align_to_sorted_bam(["cat", str(tmp_path / "in.sam")], str(bam), stderr=subprocess.DEVNULL)
        assert bam_read_count(str(bam)) == (1, 1)

    def test_failing_aligner_exits_and_removes_partial_bam(self, tmp_path):
        bam = tmp_path / "out.bam"
        with pytest.raises(SystemExit) as exc:
            align_to_sorted_bam(["false"], str(bam), stderr=subprocess.DEVNULL)
        assert exc.value.code == 1
        assert not bam.exists()
