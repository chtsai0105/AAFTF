"""Unit tests for AAFTF/assemble.py.

Covers:
  - CLI parser defaults and flag presence for the 'assemble' subcommand
  - run() guard: no left reads → sys.exit(1) for each assembler
  - spades command construction for PE, SE, and merged reads
  - spades success path: scaffolds.fasta copied to outfile
  - spades graceful handling of missing scaffolds.fasta
  - megahit command construction for PE and SE reads
  - megahit success path: final.contigs.fa copied to outfile
  - unicycler command construction
  - unimplemented methods produce a status message without crashing
"""

import sys
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.main import main

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Helper: parse 'AAFTF assemble ...' without executing the tool
# ---------------------------------------------------------------------------


def _parse_assemble(argv):
    captured = {}

    def _capture(**kwargs):
        captured["args"] = Namespace(**kwargs)

    with patch.object(sys, "argv", argv):
        with patch("aaftf.assemble.run", side_effect=_capture):
            main()
    return captured.get("args")


# ---------------------------------------------------------------------------
# Helper: build a minimal Namespace for assemble.run()
# ---------------------------------------------------------------------------


_UNSET = object()


def _make_asm_args(tmp_path, method="spades", left=_UNSET, right=None, **overrides):
    if left is _UNSET:
        left = str(tmp_path / "sample_filtered_1.fastq.gz")
    defaults = dict(
        method=method,
        left=left,
        right=right,
        longreads=None,
        merged=None,
        # Confine assemble.py's workdir (which it otherwise creates relative
        # to the CWD via a uuid/pid-based fallback name) inside tmp_path so
        # tests never leave spades_*/megahit_*/unicycler_* litter in the repo.
        workdir=str(tmp_path / f"{method}_workdir"),
        cpus=1,
        memory="16",
        careful=True,
        isolate=True,
        assembler_args=None,
        tmpdir=None,
        out=str(tmp_path / f"sample.{method}.fasta"),
        debug=False,
        pipe=True,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


# ---------------------------------------------------------------------------
# Parser defaults and required flags
# ---------------------------------------------------------------------------


class TestAssembleParser:
    def test_debug_false_by_default(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa"])
        assert args.debug is False

    def test_debug_flag_sets_true(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa", "-v"])
        assert args.debug is True

    def test_pipe_false_by_default(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa"])
        assert args.pipe is False

    def test_pipe_flag_sets_true(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa", "--pipe"])
        assert args.pipe is True

    def test_default_method_is_spades(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa"])
        assert args.method == "spades"

    def test_method_megahit(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa", "--method", "megahit"])
        assert args.method == "megahit"

    def test_method_unicycler(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa", "--method", "unicycler"])
        assert args.method == "unicycler"

    def test_default_memory_string(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa"])
        assert args.memory == "32"

    def test_custom_memory(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa", "-m", "64"])
        assert args.memory == "64"

    def test_default_cpus(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa"])
        assert args.cpus == 1

    def test_custom_cpus(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa", "-c", "8"])
        assert args.cpus == 8

    def test_careful_true_by_default(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa"])
        assert args.careful is True

    def test_no_careful_flag(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa", "--no-careful"])
        assert args.careful is False

    def test_isolate_true_by_default(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa"])
        assert args.isolate is True

    def test_no_isolate_flag(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa", "--no-isolate"])
        assert args.isolate is False

    def test_parses_left_reads(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa"])
        assert args.left == "R1.fq"

    def test_parses_right_reads(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa", "-r", "R2.fq"])
        assert args.right == "R2.fq"

    def test_right_none_by_default(self):
        args = _parse_assemble(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa"])
        assert args.right is None

    def test_assemble_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "assemble", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0


# ---------------------------------------------------------------------------
# run() guard: no left reads
# ---------------------------------------------------------------------------


class TestAssembleRunGuards:
    def test_spades_no_left_raises(self, tmp_path):
        args = _make_asm_args(tmp_path, method="spades", left=None)
        from aaftf.assemble import run

        with pytest.raises(ValueError):
            run(**vars(args))

    def test_megahit_no_left_raises(self, tmp_path):
        args = _make_asm_args(tmp_path, method="megahit", left=None)
        from aaftf.assemble import run

        with pytest.raises(ValueError):
            run(**vars(args))

    def test_unicycler_no_left_raises(self, tmp_path):
        args = _make_asm_args(tmp_path, method="unicycler", left=None)
        from aaftf.assemble import run

        with pytest.raises(ValueError):
            run(**vars(args))


# ---------------------------------------------------------------------------
# SPAdes command construction and output handling
# ---------------------------------------------------------------------------


def _run_spades(tmp_path, left, right=None, create_output=True, **extra):
    """Invoke run_spades() with subprocess mocked; return captured commands."""
    args = _make_asm_args(tmp_path, method="spades", left=left, right=right, **extra)
    cmds = []

    def _fake_run(cmd, **kw):
        cmds.append(cmd)
        if create_output and args.workdir:
            workdir = Path(args.workdir)
            workdir.mkdir(parents=True, exist_ok=True)
            (workdir / "scaffolds.fasta").write_text(">contig1\nATCGATCG\n")

    from aaftf.assemble import run_spades

    with patch("aaftf.utility.subprocess.run", side_effect=_fake_run):
        run_spades(**vars(args))
    return cmds, args


class TestAssembleRunSpades:
    def test_pe_command_includes_pe1_1(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        cmds, _ = _run_spades(tmp_path, left, right)
        assert any("--pe1-1" in c for c in cmds[0])

    def test_pe_command_includes_pe1_2(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        cmds, _ = _run_spades(tmp_path, left, right)
        assert any("--pe1-2" in c for c in cmds[0])

    def test_se_command_includes_s1(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        cmds, _ = _run_spades(tmp_path, left, right=None)
        assert any("--s1" in c for c in cmds[0])

    def test_merged_adds_s1_for_pe(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        merged = str(tmp_path / "merged.fastq.gz")
        cmds, _ = _run_spades(tmp_path, left, right, merged=merged)
        cmd_str = " ".join(cmds[0])
        assert "--s1" in cmd_str

    def test_command_includes_threads(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        cmds, _ = _run_spades(tmp_path, left, right)
        assert "--threads" in cmds[0]

    def test_command_includes_memory(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        cmds, args = _run_spades(tmp_path, left, right)
        assert "--mem" in cmds[0]
        assert args.memory in cmds[0]

    def test_careful_flag_added(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        cmds, _ = _run_spades(tmp_path, left, right, careful=True, isolate=False)
        assert "--careful" in cmds[0]

    def test_isolate_overrides_careful(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        cmds, _ = _run_spades(tmp_path, left, right, careful=False, isolate=True)
        assert "--isolate" in cmds[0]
        assert "--careful" not in cmds[0]

    def test_scaffolds_copied_to_outfile(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        out = str(tmp_path / "assembly.fasta")
        _run_spades(tmp_path, left, right, out=out, create_output=True)
        assert Path(out).exists()

    def test_missing_scaffolds_raises(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        with pytest.raises(RuntimeError, match="Spades assembly output"):
            _run_spades(tmp_path, left, right, create_output=False)

    def test_existing_workdir_restarts_from_last(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        args = _make_asm_args(tmp_path, method="spades", left=left)
        Path(args.workdir).mkdir(parents=True, exist_ok=True)
        cmds, _ = _run_spades(tmp_path, left)
        assert cmds[0][-2:] == ["--restart-from", "last"]

    def test_cov_cutoff_auto_in_command(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        cmds, _ = _run_spades(tmp_path, left, right)
        cmd_str = " ".join(cmds[0])
        assert "--cov-cutoff" in cmd_str
        assert "auto" in cmd_str


# ---------------------------------------------------------------------------
# megahit command construction and output handling
# ---------------------------------------------------------------------------


def _run_megahit(tmp_path, left, right=None, create_output=True, **extra):
    args = _make_asm_args(tmp_path, method="megahit", left=left, right=right, **extra)
    cmds = []

    def _fake_run(cmd, **kw):
        cmds.append(cmd)
        if create_output and args.workdir:
            workdir = Path(args.workdir)
            workdir.mkdir(parents=True, exist_ok=True)
            (workdir / "final.contigs.fa").write_text(">contig1\nATCGATCG\n")

    from aaftf.assemble import run_megahit

    with patch("aaftf.utility.subprocess.run", side_effect=_fake_run):
        run_megahit(**vars(args))
    return cmds, args


class TestAssembleRunMegahit:
    def test_pe_command_uses_1_2_flags(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        cmds, _ = _run_megahit(tmp_path, left, right)
        assert "-1" in cmds[0]
        assert "-2" in cmds[0]

    def test_se_command_uses_r_flag(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        cmds, _ = _run_megahit(tmp_path, left, right=None)
        assert "-r" in cmds[0]

    def test_command_includes_threads(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        cmds, _ = _run_megahit(tmp_path, left, right)
        assert "-t" in cmds[0]

    def test_final_contigs_copied_to_outfile(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        out = str(tmp_path / "assembly.megahit.fasta")
        _run_megahit(tmp_path, left, right, out=out, create_output=True)
        assert Path(out).exists()

    def test_missing_final_contigs_raises(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        right = str(tmp_path / "filtered_2.fastq.gz")
        with pytest.raises(RuntimeError, match="Megahit assembly output"):
            _run_megahit(tmp_path, left, right, create_output=False)

    def test_command_starts_with_megahit(self, tmp_path):
        left = str(tmp_path / "filtered_1.fastq.gz")
        cmds, _ = _run_megahit(tmp_path, left)
        assert cmds[0][0] == "megahit"


# ---------------------------------------------------------------------------
# unicycler command construction
# ---------------------------------------------------------------------------


def _run_unicycler(tmp_path, left, right=None, **extra):
    """Invoke run_unicycler() with subprocess mocked (writing assembly.fasta); return the command."""
    args = _make_asm_args(tmp_path, method="unicycler", left=left, right=right, **extra)
    cmds = []

    def _fake_run(cmd, **kw):
        cmds.append(cmd)
        workdir = Path(cmd[cmd.index("-o") + 1])
        workdir.mkdir(parents=True, exist_ok=True)
        (workdir / "assembly.fasta").write_text(">contig1\nATCGATCG\n")

    from aaftf.assemble import run_unicycler

    with patch("aaftf.utility.subprocess.run", side_effect=_fake_run):
        run_unicycler(**vars(args))
    return cmds[0]


class TestAssembleRunUnicycler:
    def test_pe_command_uses_short1_short2(self, tmp_path):
        cmd = _run_unicycler(tmp_path, str(tmp_path / "filtered_1.fastq.gz"), str(tmp_path / "filtered_2.fastq.gz"))
        assert "--short1" in cmd and "--short2" in cmd

    def test_pe_with_merged_keeps_pairs(self, tmp_path):
        merged = str(tmp_path / "merged.fastq.gz")
        cmd = _run_unicycler(tmp_path, str(tmp_path / "filtered_1.fastq.gz"), str(tmp_path / "filtered_2.fastq.gz"), merged=merged)
        assert "--short1" in cmd and "--short2" in cmd
        assert cmd[cmd.index("--unpaired") + 1] == merged

    def test_se_command_uses_unpaired(self, tmp_path):
        assert "--unpaired" in _run_unicycler(tmp_path, str(tmp_path / "filtered_1.fastq.gz"))

    def test_command_starts_with_unicycler(self, tmp_path):
        assert _run_unicycler(tmp_path, str(tmp_path / "filtered_1.fastq.gz"))[0] == "unicycler"


# ---------------------------------------------------------------------------
# Unimplemented / unknown methods
# ---------------------------------------------------------------------------


class TestAssembleUnimplementedMethods:
    def test_masurca_no_crash(self, tmp_path):
        args = _make_asm_args(tmp_path, method="masurca")
        from aaftf.assemble import run

        with patch("aaftf.utility.subprocess.run"):
            run(**vars(args))  # should not raise

    def test_nextdenovo_no_crash(self, tmp_path):
        args = _make_asm_args(tmp_path, method="nextdenovo")
        from aaftf.assemble import run

        with patch("aaftf.utility.subprocess.run"):
            run(**vars(args))

    def test_unknown_method_raises(self, tmp_path):
        args = _make_asm_args(tmp_path, method="unknownasm")
        from aaftf.assemble import run

        with pytest.raises(ValueError, match="Unknown assembler method"):
            run(**vars(args))
