"""Unit tests for AAFTF/mito.py helpers."""

from unittest.mock import patch

import pytest

from aaftf.mito import _rev_comp

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
