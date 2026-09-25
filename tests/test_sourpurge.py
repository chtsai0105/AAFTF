"""Unit tests for aaftf/sourpurge.py (sourmash/bwa/samtools are mocked)."""

import logging
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.sourpurge import run

pytestmark = pytest.mark.unit

_HEADER = "ID,status,superkingdom,phylum,class,order,family,genus,species,strain\n"
_FUNGUS = "Eukaryota,Ascomycota,Sordariomycetes,,,,,\n"
_BACT = "Bacteria,Proteobacteria,Gammaproteobacteria,,,,,\n"
_NOMATCH = ",,,,,,,\n"


def _setup(tmp_path, seqs):
    asm = tmp_path / "asm.fasta"
    asm.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    db = tmp_path / "lca.json"
    db.write_text("{}")
    return asm, db


def _ids(path):
    return [line[1:].strip() for line in Path(path).read_text().splitlines() if line.startswith(">")]


def _run(tmp_path, classify, bedcov=None, seqs=None, **kw):
    """Run sourpurge with the given classify CSV rows; return (kept ids, mocks)."""
    seqs = seqs or {"fungal": "A" * 100, "bact": "C" * 100, "unknown": "G" * 100}
    asm, db = _setup(tmp_path, seqs)
    out = tmp_path / "out.fasta"

    def fake_execute(cmd, cwd=None, debug=False, quiet=False):
        if cmd[:3] == ["sourmash", "lca", "classify"]:
            yield from [_HEADER, *classify]
        else:
            yield from bedcov or []

    with (
        patch("aaftf.sourpurge.run_cmd") as rc,
        patch("aaftf.sourpurge.execute", side_effect=fake_execute) as ex,
        patch("aaftf.sourpurge.align_to_sorted_bam") as bam,
    ):
        kw.setdefault("phylum", ["Ascomycota"])
        run(input=str(asm), outfile=str(out), sourdb=str(db), workdir=str(tmp_path / "wd"), **kw)
    return (_ids(out) if out.exists() else None), rc, ex, bam


def test_taxonomy_filter_keeps_phylum_and_nomatch(tmp_path):
    rows = ["fungal,found," + _FUNGUS, "bact,found," + _BACT, "unknown,nomatch," + _NOMATCH]
    kept, rc, ex, bam = _run(tmp_path, rows)
    assert kept == ["fungal", "unknown"]
    assert rc.call_args_list[0].args[0][:2] == ["sourmash", "compute"]
    assert ex.call_args_list[0].args[0][:5] == ["sourmash", "lca", "classify", "--db", str(tmp_path / "lca.json")]
    rc.assert_called_once()  # no bwa index without reads
    bam.assert_not_called()
    assert (tmp_path / "asm.sourmash-taxonomy.csv").read_text().startswith(_HEADER)


def test_status_must_be_exactly_found(tmp_path):
    """Only status "found" counts as a classification; a status merely containing it is ignored."""
    kept, *_ = _run(tmp_path, ["fungal,found," + _FUNGUS, "bact,notfound," + _BACT])
    assert kept == ["fungal", "bact", "unknown"]


def test_phylum_matches_any_rank_exactly(tmp_path):
    # matching is exact membership in the lineage list, so a class name works and a prefix does not
    rows = ["fungal,found," + _FUNGUS, "bact,found," + _BACT]
    kept, *_ = _run(tmp_path, rows, phylum=["Gammaproteobacteria", "Ascomy"])
    assert kept == ["bact", "unknown"]


def test_taxonomy_only_returns_early(tmp_path):
    kept, rc, *_ = _run(tmp_path, ["bact,found," + _BACT], taxonomy=True)
    assert kept is None
    assert not (tmp_path / "asm.sourmash-taxonomy.csv").exists()


def _bedcov(covs):
    return [f"{k}\t0\t{length}\t{int(length * c)}\n" for k, (length, c) in covs.items()]


@pytest.mark.parametrize("paired", [False, True])
def test_coverage_filter_with_reads(tmp_path, paired):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    seqs = {"big": "A" * 1000, "mid": "C" * 500, "low": "G" * 400, "edge": "T" * 300}
    # N50 = 500 bp, so N50 contigs are big+mid (avg 100X): cutoff 5% -> 5X; <= cutoff dropped
    bed = _bedcov({"big": (1000, 100), "mid": (500, 100), "low": (400, 5.5), "edge": (300, 5)})
    kw = {"read1": str(r1)}
    if paired:
        kw["read2"] = str(r2)
    kept, rc, ex, bam = _run(tmp_path, [], bedcov=bed, seqs=seqs, **kw)
    assert kept == ["big", "mid", "low"]
    assert rc.call_args_list[1].args[0] == ["bwa", "index", "sourmashed-tax-screen.fasta"]
    bwa = bam.call_args.args[0]
    assert bwa[:5] == ["bwa", "mem", "-t", "1", "sourmashed-tax-screen.fasta"]
    assert bwa[5:] == ([str(r1.resolve()), str(r2.resolve())] if paired else [str(r1.resolve())])
    assert ex.call_args_list[1].args[0] == ["samtools", "bedcov", "assembly.bed", "remapped.bam"]


def test_mincovpct_changes_cutoff(tmp_path):
    seqs = {"big": "A" * 1000, "mid": "C" * 500}
    kept, *_ = _run(tmp_path, [], bedcov=_bedcov({"big": (1000, 100), "mid": (500, 20)}), seqs=seqs, read1=str(tmp_path / "r1.fq"), mincovpct=20)
    assert kept == ["big"]


def test_no_n50_coverage_warns_and_skips(tmp_path, caplog):
    caplog.set_level(logging.INFO)
    seqs = {"big": "A" * 1000, "low": "G" * 400}
    kept, *_ = _run(tmp_path, [], bedcov=_bedcov({"low": (400, 0)}), seqs=seqs, read1=str(tmp_path / "r1.fq"))
    assert kept == ["big", "low"]
    assert "skipping the low-coverage filter" in caplog.text


def test_existing_bam_skips_alignment(tmp_path):
    (tmp_path / "wd").mkdir()
    (tmp_path / "wd" / "remapped.bam").write_text("")
    kept, rc, ex, bam = _run(tmp_path, [], bedcov=_bedcov({"fungal": (100, 10)}), seqs={"fungal": "A" * 100}, read1=str(tmp_path / "r1.fq"))
    bam.assert_not_called()
    assert kept == ["fungal"]


def test_default_db_from_sourdb_type(tmp_path):
    asm, _ = _setup(tmp_path, {"a": "A" * 10})
    with (
        patch("aaftf.sourpurge.require_databases", return_value=[str(tmp_path / "gtdb.json")]) as req,
        patch("aaftf.sourpurge.run_cmd"),
        patch("aaftf.sourpurge.execute", return_value=iter([])) as ex,
    ):
        run(input=str(asm), outfile=str(tmp_path / "o.fa"), phylum=["X"], sourdb_type="gtdbrep", workdir=str(tmp_path / "wd"))
    assert req.call_args.args[0] == ["sm_gtdbrep"]
    assert str(tmp_path / "gtdb.json") in ex.call_args.args[0]


def test_pipe_hides_hint(tmp_path, caplog):
    caplog.set_level(logging.INFO)
    _run(tmp_path, [], pipe=True)
    assert "AAFTF rmdup" not in caplog.text
    _run(tmp_path, [], pipe=False)
    assert "AAFTF rmdup" in caplog.text
