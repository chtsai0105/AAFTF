"""Unit tests for aaftf/pipeline.py: every step gets its own CLI defaults plus the pipeline's options."""

import argparse as ap
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf import pipeline
from aaftf._menu import register_subcommands

pytestmark = pytest.mark.unit

STEPS = ["trim", "filter", "assemble", "vecscreen", "sourpurge", "rmdup", "polish", "sort", "assess"]
_MODULE = {"filter": "aaftf_filter", "sort": "aaftf_sort"}
# the kwarg that names each step's output file (trim and filter derive theirs from basename)
_OUTPUT = {"assemble": "out", "vecscreen": "outfile", "sourpurge": "outfile", "rmdup": "out", "polish": "outfile", "sort": "out"}


def _cli_defaults(step):
    """Return ``{dest: default}`` for ``AAFTF <step>``."""
    parser = register_subcommands(ap.ArgumentParser()).choices[step]
    return {a.dest: a.default for a in parser._actions if a.dest != "help"}


def _run_pipeline(tmp_path, read2=True, **options):
    """Run pipeline.run() with every step's run() replaced; return {step: kwargs it was called with}."""
    base = str(tmp_path / "sample")
    calls = {}

    def _fake(step):
        def _run(**kwargs):
            calls[step] = kwargs
            outputs = {"trim": [f"{base}_1P.fastq.gz"], "filter": [f"{base}_filtered_1.fastq.gz"]}.get(step, [kwargs.get(_OUTPUT.get(step, ""))])
            for out in filter(None, outputs):
                Path(out).write_text("data")

        return _run

    patches = [patch(f"aaftf.pipeline.{_MODULE.get(step, step)}.run", side_effect=_fake(step)) for step in STEPS]
    for p in patches:
        p.start()
    try:
        pipeline.run(read1="R1.fq.gz", read2="R2.fq.gz" if read2 else None, basename=base, phylum=["Ascomycota"], **options)
    finally:
        for p in patches:
            p.stop()
    return calls


class TestSteps:
    def test_runs_every_step_and_no_optional_ones(self, tmp_path):
        assert list(_run_pipeline(tmp_path)) == STEPS

    @pytest.mark.parametrize("step", STEPS)
    def test_step_gets_every_cli_option(self, tmp_path, step):
        assert set(_run_pipeline(tmp_path)[step]) == set(_cli_defaults(step)) | {"pipe"}

    @pytest.mark.parametrize("step", STEPS)
    def test_step_hints_suppressed(self, tmp_path, step):
        assert _run_pipeline(tmp_path)[step]["pipe"] is True

    def test_existing_output_skips_step(self, tmp_path):
        Path(f"{tmp_path / 'sample'}_1P.fastq.gz").write_text("done")
        assert "trim" not in _run_pipeline(tmp_path)

    def test_missing_output_raises(self, tmp_path):
        with patch("aaftf.pipeline.trim.run"):
            with pytest.raises(RuntimeError, match="AAFTF trim failed"):
                pipeline.run(read1="R1.fq.gz", basename=str(tmp_path / "sample"), phylum=["Ascomycota"])


class TestDefaults:
    """Options the pipeline does not expose keep each step's own default."""

    @pytest.mark.parametrize(
        "step, option",
        [
            ("trim", "method"),
            ("trim", "avgqual"),
            ("filter", "aligner"),
            ("filter", "screen_local"),
            ("assemble", "isolate"),
            ("assemble", "careful"),
            ("assemble", "merged"),
            ("vecscreen", "percent_id"),
            ("vecscreen", "stringency"),
            ("sourpurge", "kmer"),
            ("sourpurge", "sourdb_type"),
            ("rmdup", "percent_id"),
            ("rmdup", "percent_cov"),
            ("polish", "method"),
            ("sort", "name"),
            ("assess", "telomere_monomer"),
            ("assess", "report"),
        ],
    )
    def test_step_default_kept(self, tmp_path, step, option):
        assert _run_pipeline(tmp_path)[step][option] == _cli_defaults(step)[option]

    @pytest.mark.parametrize("step", ["trim", "filter", "assemble", "polish"])
    def test_memory_default_kept_without_pipeline_memory(self, tmp_path, step):
        assert _run_pipeline(tmp_path)[step]["memory"] == _cli_defaults(step)["memory"]


class TestPipelineOptions:
    """Options given to the pipeline override the step defaults."""

    def test_memory_passed_to_each_step(self, tmp_path):
        calls = _run_pipeline(tmp_path, memory=12)
        assert [calls[s]["memory"] for s in ("trim", "filter", "assemble", "polish")] == [12, 12, 12, 12]

    @pytest.mark.parametrize("step", ["trim", "filter", "assemble", "vecscreen", "sourpurge", "rmdup", "polish"])
    def test_cpus_passed(self, tmp_path, step):
        assert _run_pipeline(tmp_path, cpus=6)[step]["cpus"] == 6

    @pytest.mark.parametrize("step", ["filter", "assemble", "vecscreen", "sourpurge", "rmdup", "polish"])
    def test_workdir_passed(self, tmp_path, step):
        assert _run_pipeline(tmp_path, workdir="wd")[step]["workdir"] == "wd"

    @pytest.mark.parametrize("step", STEPS)
    def test_debug_passed(self, tmp_path, step):
        assert _run_pipeline(tmp_path, debug=True)[step]["debug"] is True

    def test_step_specific_options(self, tmp_path):
        calls = _run_pipeline(tmp_path, minlen=50, mincontiglen=1000, method="megahit", mincovpct=7, sourdb="db.zip", screen_accessions=["NC_1"])
        assert calls["trim"]["minlen"] == 50
        assert calls["rmdup"]["minlen"] == calls["sort"]["minlen"] == 1000
        assert calls["assemble"]["method"] == "megahit"
        assert (calls["sourpurge"]["mincovpct"], calls["sourpurge"]["sourdb"]) == (7, "db.zip")
        assert calls["filter"]["screen_accessions"] == ["NC_1"]

    def test_pipeline_method_is_only_the_assembler(self, tmp_path):
        calls = _run_pipeline(tmp_path, method="megahit")
        assert calls["trim"]["method"] == _cli_defaults("trim")["method"]
        assert calls["polish"]["method"] == _cli_defaults("polish")["method"]


class TestFileChaining:
    def test_each_step_reads_the_previous_output(self, tmp_path):
        calls = _run_pipeline(tmp_path)
        base = str(tmp_path / "sample")
        assert (calls["filter"]["read1"], calls["filter"]["read2"]) == (f"{base}_1P.fastq.gz", f"{base}_2P.fastq.gz")
        assert calls["assemble"]["read1"] == f"{base}_filtered_1.fastq.gz"
        assert calls["vecscreen"]["infile"] == f"{base}.spades.fasta"
        assert calls["sourpurge"]["input"] == f"{base}.vecscreen.fasta"
        assert calls["rmdup"]["input"] == f"{base}.sourpurge.fasta"
        assert calls["polish"]["infile"] == f"{base}.rmdup.fasta"
        assert calls["sort"]["input"] == f"{base}.polish.fasta"
        assert calls["assess"]["input"] == f"{base}.final.fasta"

    def test_single_end_reads_have_no_read2(self, tmp_path):
        calls = _run_pipeline(tmp_path, read2=False)
        assert all(calls[s]["read2"] is None for s in ("trim", "filter", "assemble", "sourpurge", "polish"))
