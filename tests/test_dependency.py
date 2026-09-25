"""Unit tests for aaftf/dependency.py."""

import pytest

from aaftf.dependency import _print_tool_table

pytestmark = pytest.mark.unit

RESULTS = [("bwa", "/opt/bin/bwa", "filter, depth"), ("bowtie2", None, "filter")]


def test_full_table_shows_location_and_users(capsys):
    assert _print_tool_table("Checking...", RESULTS) == ["bowtie2"]
    out = capsys.readouterr().out
    assert "/opt/bin/bwa" in out and "required by: filter" in out


def test_quiet_table_shows_only_status_and_name(capsys):
    assert _print_tool_table("Checking...", RESULTS, quiet=True) == ["bowtie2"]
    assert capsys.readouterr().out.splitlines() == ["  [OK]      bwa", "  [MISSING] bowtie2"]
