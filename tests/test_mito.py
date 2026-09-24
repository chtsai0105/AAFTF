"""Unit tests for AAFTF/mito.py helpers."""

from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.mito import _rev_comp, run

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# _rev_comp
# ---------------------------------------------------------------------------


class TestRevComp:
    def test_simple(self):
        assert _rev_comp("ATCG") == "CGAT"

    def test_complement_only(self):
        assert _rev_comp("AAAA") == "TTTT"
        assert _rev_comp("CCCC") == "GGGG"

    def test_palindrome(self):
        assert _rev_comp("AATTAATT") == "AATTAATT"

    def test_preserves_case(self):
        assert _rev_comp("atcg") == "cgat"
        assert _rev_comp("AAcgTT") == "AAcgTT"

    def test_iupac_codes(self):
        assert _rev_comp("RYKMN") == "NKMRY"

    def test_longer_sequence(self):
        assert _rev_comp("ATCGATCG") == "CGATCGAT"

    def test_all_bases(self):
        # A↔T, C↔G
        assert _rev_comp("ACGT") == "ACGT"


class TestOrientToStart:
    """_orient_to_start rotates the circular genome so the start gene's alignment begins at position 0."""

    def _orient(self, tmp_path, hits, seq):
        from aaftf.mito import _orient_to_start
        from aaftf.utility import PafHit

        fasta_in, fasta_out = tmp_path / "in.fa", tmp_path / "out.fa"
        fasta_in.write_text(f">mt\n{seq}\n")
        paf = [PafHit("COB", 100, qs, 100, strand, "mt", len(seq), ts, te, 100, 100, 60) for qs, strand, ts, te in hits]
        with patch("aaftf.mito.paf_hits", return_value=iter(paf)):
            _orient_to_start(str(fasta_in), str(fasta_out), folder=str(tmp_path))
        return "".join(fasta_out.read_text().splitlines()[1:])

    def test_forward_hit_rotates_to_target_start(self, tmp_path):
        seq = "AAAACCCCGGGGTTTT"
        assert self._orient(tmp_path, [(0, "+", 4, 8)], seq) == "CCCCGGGGTTTTAAAA"

    def test_reverse_hit_rotates_and_reverse_complements(self, tmp_path):
        seq = "AAAACCCCGGGGTTTT"
        rotated = seq[8:] + seq[:8]
        assert self._orient(tmp_path, [(0, "-", 4, 8)], seq) == _rev_comp(rotated)

    def test_no_hit_leaves_sequence_unrotated(self, tmp_path):
        seq = "AAAACCCCGGGGTTTT"
        assert self._orient(tmp_path, [], seq) == seq


class TestNovoplastyInputs:
    """run() writes the NOVOPlasty config from the bundled template and the bundled seed (package data)."""

    def _run(self, tmp_path, **kwargs):
        from aaftf.mito import run

        workdir = tmp_path / "wd"
        with patch("aaftf.mito.require_tools"), patch("aaftf.mito.estimate_read_length", return_value=150):
            with patch("aaftf.mito.subprocess.Popen") as popen:
                popen.return_value.communicate.return_value = (b"", b"")
                # no NOVOPlasty output is faked, so run() stops after writing its inputs
                with pytest.raises(RuntimeError, match="did not produce an assembly"):
                    run(left=str(tmp_path / "R1.fq"), right=str(tmp_path / "R2.fq"), out=str(tmp_path / "mt.fa"), workdir=str(workdir), **kwargs)
        return workdir

    def test_default_seed_copied_from_package_and_config_filled(self, tmp_path):
        workdir = self._run(tmp_path, memory=6)
        seed = workdir / "mito-seed.fasta"
        assert seed.read_text().startswith(">")
        config = (workdir / "novo-config.txt").read_text()
        assert "<" not in config.replace("<=", "")  # every <PLACEHOLDER> replaced
        assert str(seed.resolve()) in config and "150" in config and "6" in config

    def test_explicit_seed_not_copied(self, tmp_path):
        own_seed = tmp_path / "myseed.fa"
        own_seed.write_text(">s\nACGT\n")
        workdir = self._run(tmp_path, seed=str(own_seed))
        assert not (workdir / "mito-seed.fasta").exists()
        assert str(own_seed.resolve()) in (workdir / "novo-config.txt").read_text()


class TestNovoplastyOutputChoice:
    """run() picks NOVOPlasty's circular assembly first, then Contigs_1, then the uncircularised file."""

    def _run(self, tmp_path, outputs):
        def _fake_popen(cmd, cwd, **kwargs):
            for name, seq in outputs.items():
                Path(cwd, name).write_text(f">x\n{seq}\n")
            return type("P", (), {"communicate": lambda self: None})()

        orient = []
        out = tmp_path / "mt.fasta"
        with patch("aaftf.mito.require_tools"), patch("aaftf.mito.estimate_read_length", return_value=150):
            with patch("aaftf.mito.subprocess.Popen", side_effect=_fake_popen):
                with patch("aaftf.mito._orient_to_start", side_effect=lambda draft, out, **kw: orient.append(Path(draft).name)):
                    run(left="R1.fq", right="R2.fq", out=str(out), workdir=str(tmp_path / "wd"))
        return orient, out

    def test_circular_preferred(self, tmp_path):
        orient, _ = self._run(tmp_path, {"Contigs_1_x.fasta": "AAAA", "Circularized_assembly_1_x.fasta": "CCCC", "Uncircularized_assemblies_x.fasta": "GGGG"})
        assert orient == ["Circularized_assembly_1_x.fasta"]

    def test_contigs_preferred_over_uncircularized(self, tmp_path):
        orient, out = self._run(tmp_path, {"Uncircularized_assemblies_x.fasta": "GGGG", "Contigs_1_x.fasta": "AAAA"})
        assert orient == [] and "AAAA" in out.read_text()

    def test_no_assembly_raises(self, tmp_path):
        with pytest.raises(RuntimeError, match="NOVOplasty did not produce an assembly"):
            self._run(tmp_path, {})
