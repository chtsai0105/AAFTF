"""Unit tests for aaftf/fcs_gx_purge.py (run_gx.py is mocked)."""

import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.fcs_gx_purge import run

pytestmark = pytest.mark.unit

_HEADER = "##[header]\n#seq_id\tstart_pos\tend_pos\tseq_len\taction\tdiv\tagg_cont_cov\ttop_tax_name\n"


def _run(tmp_path, report_rows, returncode=0, write_report=True):
    """Run fcs_gx_purge.run() on a 3-contig assembly with a fake run_gx.py; return the kept contig IDs."""
    asm = tmp_path / "asm.fasta"
    asm.write_text(">keep\nACGT\n>bad\nACGT\n>trimme\nACGT\n")
    db = tmp_path / "gxdb" / "all"
    db.parent.mkdir()
    Path(f"{db}.gxi").write_text("")
    workdir = tmp_path / "wd"
    out = tmp_path / "out.fasta"

    def _fake_run(cmd, **kwargs):
        if write_report:
            Path(workdir, "asm.4890.fcs_gx_report.txt").write_text(_HEADER + "".join("\t".join(r) + "\n" for r in report_rows))
        return subprocess.CompletedProcess(cmd, returncode)

    with patch("aaftf.fcs_gx_purge.subprocess.run", side_effect=_fake_run):
        run(input=str(asm), outfile=str(out), workdir=str(workdir), db=str(db), pipe=True)
    return [line[1:] for line in out.read_text().split() if line.startswith(">")]


def test_only_exclude_rows_are_dropped(tmp_path, caplog):
    rows = [["bad", "1", "4", "4", "EXCLUDE", "anml:insects", "100", "fly"], ["trimme", "1", "2", "4", "TRIM", "prok:bact", "50", "E. coli"]]
    assert _run(tmp_path, rows) == ["keep", "trimme"]
    assert "TRIM for trimme:1..2" in caplog.text


def test_action_prefix_accepted(tmp_path):
    assert _run(tmp_path, [["bad", "1", "4", "4", "ACTION_EXCLUDE", "x", "100", "y"]]) == ["keep", "trimme"]


def test_failed_run_gx_raises(tmp_path):
    with pytest.raises(RuntimeError, match="run_gx.py failed"):
        _run(tmp_path, [], returncode=1)


def test_missing_report_raises(tmp_path):
    with pytest.raises(RuntimeError, match="did not write"):
        _run(tmp_path, [], write_report=False)
