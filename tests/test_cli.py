"""Tests for the AAFTF CLI framework.

Covers subcommand registration, argparse defaults, and required-argument
enforcement. Each subcommand parser binds its own subtool's run() function
directly via parser.set_defaults(func=<module>.run) in AAFTF/_menu.py, and
main() invokes it as args.func(**vars(args)); these tests intercept that by
patching the relevant module's run() before it fires.

No external bioinformatics tools are invoked.
"""

import sys
from argparse import Namespace
from unittest.mock import patch

import pytest

from aaftf.main import main

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Subcommand help (argparse –– no tool execution)
# ---------------------------------------------------------------------------

# Maps the subcommand used in this file's test argv to the module whose
# run() must be intercepted before it actually executes.
_COMMAND_MODULE = {
    "depth": "aaftf.depth",
    "assess": "aaftf.assess",
    "sort": "aaftf.sort",
    "fix_tbl": "aaftf.fix_tbl",
    "polish": "aaftf.polish",
    "assemble": "aaftf.assemble",
}


def _parse_with_main(argv):
    """Run main() with sys.argv patched, intercept the subtool's run() before it fires."""
    captured = {}

    def _capture(**kwargs):
        captured["args"] = Namespace(**kwargs)

    module = _COMMAND_MODULE[argv[1]]
    with patch.object(sys, "argv", argv):
        with patch(f"{module}.run", side_effect=_capture):
            main()
    return captured.get("args")


class TestSubcommandRegistration:
    def test_top_level_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0

    def test_depth_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "depth", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0

    def test_assess_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "assess", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0

    def test_sort_help_exits_zero(self):
        with patch.object(sys, "argv", ["AAFTF", "sort", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0


# ---------------------------------------------------------------------------
# depth parser — argument parsing and defaults
# ---------------------------------------------------------------------------


class TestDepthParser:
    def test_parses_input_flag(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "genome.fa", "-l", "left.fq"])
        assert args.input == "genome.fa"

    def test_parses_left_reads(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "left.fq"])
        assert args.left == "left.fq"

    def test_parses_right_reads(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-r", "r.fq"])
        assert args.right == "r.fq"

    def test_parses_longreads(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-lr", "lr.fq"])
        assert args.longreads == "lr.fq"

    def test_default_out(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.out == "coverage_stats.txt"

    def test_custom_out(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-o", "my_report.txt"])
        assert args.out == "my_report.txt"

    def test_default_cpus(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.cpus == 1

    def test_custom_cpus(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-c", "8"])
        assert args.cpus == 8

    def test_default_aligner(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.aligner == "minimap2"

    def test_bwa_aligner(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "--aligner", "bwa"])
        assert args.aligner == "bwa"

    def test_longread_preset_none_by_default(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-lr", "lr.fq"])
        assert args.longread_preset is None

    def test_longread_preset_map_pb(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-lr", "lr.fq", "--longread_preset", "map-pb"])
        assert args.longread_preset == "map-pb"

    def test_debug_flag_false_by_default(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.debug is False

    def test_debug_flag_set(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-v"])
        assert args.debug is True

    def test_pipe_flag_false_by_default(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.pipe is False

    def test_missing_input_exits_nonzero(self):
        with patch.object(sys, "argv", ["AAFTF", "depth", "-l", "l.fq"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code != 0

    def test_default_plot_format_is_pdf(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.plot_format == "pdf"

    def test_plot_format_svg(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "--plot-format", "svg"])
        assert args.plot_format == "svg"

    def test_plot_format_png(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "--plot-format", "png"])
        assert args.plot_format == "png"

    def test_no_plot_false_by_default(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq"])
        assert args.no_plot is False

    def test_no_plot_flag(self):
        args = _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "--no-plot"])
        assert args.no_plot is True


# ---------------------------------------------------------------------------
# assess parser — common flags present
# ---------------------------------------------------------------------------


class TestAssessParser:
    def test_has_debug_flag(self):
        args = _parse_with_main(["AAFTF", "assess", "-i", "g.fa", "-v"])
        assert args.debug is True

    def test_has_pipe_flag(self):
        args = _parse_with_main(["AAFTF", "assess", "-i", "g.fa", "--pipe"])
        assert args.pipe is True

    def test_default_telomere_monomer(self):
        args = _parse_with_main(["AAFTF", "assess", "-i", "g.fa"])
        assert args.telomere_monomer == "TAAC{3,5}"

    def test_default_telomere_n_repeat(self):
        args = _parse_with_main(["AAFTF", "assess", "-i", "g.fa"])
        assert args.telomere_n_repeat == 2


# ---------------------------------------------------------------------------
# sort parser — common flags present
# ---------------------------------------------------------------------------


class TestSortParser:
    def test_has_debug_flag(self):
        args = _parse_with_main(["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa", "-v"])
        assert args.debug is True

    def test_has_pipe_flag(self):
        args = _parse_with_main(["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa", "--pipe"])
        assert args.pipe is True

    def test_default_name_prefix(self):
        args = _parse_with_main(["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa"])
        assert args.name == "scaffold"


# ---------------------------------------------------------------------------
# fix_tbl parser — common flags present (added in audit 2026-05-02)
# ---------------------------------------------------------------------------


class TestFixTblParser:
    def test_has_debug_flag(self, tmp_path):
        tbl = tmp_path / "in.tbl"
        rpt = tmp_path / "rpt.csv"
        out = tmp_path / "out.tbl"
        tbl.write_text("")
        rpt.write_text("")
        args = _parse_with_main(["AAFTF", "fix_tbl", "-t", str(tbl), "-r", str(rpt), "-o", str(out), "-v"])
        assert args.debug is True

    def test_has_pipe_flag(self, tmp_path):
        tbl = tmp_path / "in.tbl"
        rpt = tmp_path / "rpt.csv"
        out = tmp_path / "out.tbl"
        tbl.write_text("")
        rpt.write_text("")
        args = _parse_with_main(["AAFTF", "fix_tbl", "-t", str(tbl), "-r", str(rpt), "-o", str(out), "--pipe"])
        assert args.pipe is True

    def test_debug_false_by_default(self, tmp_path):
        tbl = tmp_path / "in.tbl"
        rpt = tmp_path / "rpt.csv"
        out = tmp_path / "out.tbl"
        tbl.write_text("")
        rpt.write_text("")
        args = _parse_with_main(["AAFTF", "fix_tbl", "-t", str(tbl), "-r", str(rpt), "-o", str(out)])
        assert args.debug is False

    def test_pipe_false_by_default(self, tmp_path):
        tbl = tmp_path / "in.tbl"
        rpt = tmp_path / "rpt.csv"
        out = tmp_path / "out.tbl"
        tbl.write_text("")
        rpt.write_text("")
        args = _parse_with_main(["AAFTF", "fix_tbl", "-t", str(tbl), "-r", str(rpt), "-o", str(out)])
        assert args.pipe is False


# ---------------------------------------------------------------------------
# assess parser — telomere_window argument (added in audit 2026-05-02)
# ---------------------------------------------------------------------------


class TestAssessWindowParser:
    def test_default_telomere_window(self):
        args = _parse_with_main(["AAFTF", "assess", "-i", "g.fa"])
        assert args.telomere_window == 200

    def test_custom_telomere_window(self):
        args = _parse_with_main(["AAFTF", "assess", "-i", "g.fa", "--telomere_window", "500"])
        assert args.telomere_window == 500


class TestCpusCappedToAvailable:
    """main() lowers -c/--cpus to the CPUs available to the job, with a warning."""

    def _parse(self, cpus, available):
        with patch("aaftf.main.available_cpus", return_value=available):
            return _parse_with_main(["AAFTF", "depth", "-i", "g.fa", "-l", "l.fq", "-c", str(cpus)])

    # main() configures logging to stderr (propagation off), so read captured stderr, not caplog
    def test_too_many_cpus_capped(self, capsys):
        assert self._parse(cpus=8, available=2).cpus == 2
        assert "WARNING: -c/--cpus 8 is more than the 2 CPUs available" in capsys.readouterr().err

    def test_within_available_unchanged(self, capsys):
        assert self._parse(cpus=2, available=4).cpus == 2
        assert "CPUs available" not in capsys.readouterr().err


class TestMemoryCappedToAvailable:
    """main() lowers -m/--memory to the RAM available, keeping the option's type."""

    def test_too_much_memory_capped(self, capsys):
        with patch("aaftf.main.get_ram", return_value=5.6):
            args = _parse_with_main(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq", "-m", "16"])
        assert args.memory == 5
        assert "WARNING: -m/--memory 16 GB is more than the 5.6 GB of RAM available" in capsys.readouterr().err

    def test_string_memory_stays_a_string(self):
        with patch("aaftf.main.get_ram", return_value=10.0):
            args = _parse_with_main(["AAFTF", "assemble", "-l", "R1.fq", "-o", "out.fa"])
        assert args.memory == "10"

    def test_within_available_unchanged(self, capsys):
        with patch("aaftf.main.get_ram", return_value=64.0):
            args = _parse_with_main(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.memory == 16
        assert "RAM available" not in capsys.readouterr().err

    def test_never_below_1gb(self):
        with patch("aaftf.main.get_ram", return_value=0.3):
            args = _parse_with_main(["AAFTF", "polish", "-i", "asm.fa", "-l", "R1.fq"])
        assert args.memory == 1


class TestVersionAndVerboseFlags:
    @pytest.mark.parametrize("flag", ["-V", "--version"])
    def test_version_flag_prints_version(self, flag, capsys):
        with patch.object(sys, "argv", ["AAFTF", flag]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code == 0
        assert capsys.readouterr().out.startswith("AAFTF ")

    def test_verbose_long_flag_sets_debug(self):
        args = _parse_with_main(["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa", "--verbose"])
        assert args.debug is True

    def test_old_debug_flag_rejected(self):
        with patch.object(sys, "argv", ["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa", "--debug"]):
            with pytest.raises(SystemExit) as exc:
                main()
        assert exc.value.code != 0


class TestMainExitCodes:
    """main() turns subcommand exceptions into shell exit codes instead of tracebacks."""

    def _run(self, side_effect):
        with patch.object(sys, "argv", ["AAFTF", "sort", "-i", "in.fa", "-o", "out.fa"]):
            with patch("aaftf.sort.run", side_effect=side_effect):
                return main()

    def test_success_returns_0(self):
        assert self._run(None) == 0

    def test_missing_file_returns_2(self):
        assert self._run(FileNotFoundError("no such file")) == 2

    def test_other_error_returns_1(self, capsys):
        assert self._run(RuntimeError("tool failed")) == 1
        assert "tool failed" in capsys.readouterr().err

    def test_keyboard_interrupt_returns_130(self):
        assert self._run(KeyboardInterrupt()) == 130

    def test_no_arguments_returns_1(self):
        with patch.object(sys, "argv", ["AAFTF"]):
            assert main() == 1
