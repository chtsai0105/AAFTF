"""Unit tests for AAFTF/mito.py helpers."""

import pytest

from AAFTF.mito import _rev_comp

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
