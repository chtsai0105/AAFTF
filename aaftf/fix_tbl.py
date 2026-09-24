"""Fix an NCBI .tbl annotation table after NCBI FCS trimmed or excluded contigs.

Feature coordinates are shifted by the bases trimmed off the start of each
contig, features are clipped at the new contig ends (and marked partial with
``<``/``>`` where they were cut), and features that fell entirely inside a
trimmed region, or on an excluded contig, are dropped.
"""

import logging
import re
from dataclasses import dataclass, field

__all__ = ["TblFeature", "run", "fix_tbl", "parse_tbl", "parse_adjustments"]


logger = logging.getLogger(__name__)

_HEADER_RE = re.compile(r">Feature\s+(\S+)")
_INTERVAL_RE = re.compile(r"^([<>]?\d+)\t([<>]?\d+)(?:\t(\S+))?")
_RANGE_RE = re.compile(r"(\d+)\.\.(\d+)")


@dataclass
class TblFeature:
    """One feature of a .tbl file: its key, coordinate intervals and qualifier lines.

    Attributes:
        key: Feature key, e.g. ``gene`` or ``CDS``.
        intervals: ``[start, end]`` coordinate strings as written in the file,
            including any ``<``/``>`` partial markers; ``start > end`` on the minus strand.
        qualifiers: Qualifier lines exactly as read (with their leading tabs).
    """

    key: str
    intervals: list = field(default_factory=list)
    qualifiers: list = field(default_factory=list)

    def to_tbl(self):
        """Return the feature as .tbl text."""
        lines = [f"{self.intervals[0][0]}\t{self.intervals[0][1]}\t{self.key}"]
        lines += [f"{start}\t{end}" for start, end in self.intervals[1:]]
        lines += self.qualifiers
        return "\n".join(lines) + "\n"


def run(table, report, output, **kwargs):
    """Run the fix_tbl subcommand."""
    fix_tbl(table, report, output)


def fix_tbl(tbl_fh, adjustment_fh, output_handle):
    """Write ``tbl_fh`` to ``output_handle`` with coordinates fixed for the FCS trims in ``adjustment_fh``."""
    features = parse_tbl(tbl_fh)
    adjustments = parse_adjustments(adjustment_fh)
    for seqid, feats in features.items():
        window = _kept_window(seqid, adjustments.get(seqid, []))
        if window is None:
            logger.warning(f"{seqid} was excluded by FCS; dropping its {len(feats)} feature(s)")
            continue
        print(f">Feature {seqid}", file=output_handle)
        for feat in feats:
            fixed = _clip_feature(feat, *window)
            if fixed is None:
                logger.warning(f"dropping {feat.key} at {seqid}:{feat.intervals[0][0]}..{feat.intervals[-1][1]}; it lies inside a trimmed region")
                continue
            output_handle.write(fixed.to_tbl())


def parse_tbl(tbl_file_handle):
    """Parse an NCBI .tbl file into ``{seqid: [TblFeature, ...]}`` (in file order)."""
    features = {}
    sequence_name = None
    for lineno, line in enumerate(tbl_file_handle, 1):
        line = line.rstrip("\r\n")
        if not line.strip() or line.startswith("#"):
            continue
        if line.startswith(">"):
            m = _HEADER_RE.match(line)
            if not m:
                raise ValueError(f"line {lineno}: unexpected header in tbl file: {line}")
            sequence_name = m.group(1)
            features[sequence_name] = []
            continue
        if sequence_name is None:
            raise ValueError(f"line {lineno}: feature line before the first '>Feature' header: {line}")
        feats = features[sequence_name]
        if m := _INTERVAL_RE.match(line):
            start, end, key = m.groups()
            if key:
                feats.append(TblFeature(key, [[start, end]]))
            elif feats and not feats[-1].qualifiers:
                feats[-1].intervals.append([start, end])  # next interval of a multi-interval feature
            else:
                raise ValueError(f"line {lineno}: coordinate line without a feature key: {line}")
        elif line.startswith("\t") and feats:
            feats[-1].qualifiers.append(line)
        else:
            logger.warning(f"line {lineno}: skipping unrecognised tbl line: {line}")
    return features


def parse_adjustments(adj_file_handle):
    """Parse an NCBI FCS action report into ``{seqid: [(length, action, start, end), ...]}``.

    The report is tab-separated (accession, length, action, range(s), ...) after a
    ``#accession`` header line; lines before that header are ignored. An
    ``ACTION_EXCLUDE`` row without ranges covers the whole sequence.
    """
    adjustments = {}
    in_table = False
    for line in adj_file_handle:
        if not in_table:
            in_table = line.startswith("#accession")
            continue
        row = line.rstrip("\r\n").split("\t")
        if len(row) < 3 or not row[0] or row[0].startswith("#"):
            continue
        seqid, length, action = row[0], row[1], row[2]
        try:
            length = int(length)
        except ValueError:
            logger.warning(f"skipping FCS report line with a non-numeric length: {line.strip()}")
            continue
        ranges = [(int(s), int(e)) for s, e in _RANGE_RE.findall(row[3] if len(row) > 3 else "")]
        if not ranges and action == "ACTION_EXCLUDE":
            ranges = [(1, length)]
        adjustments.setdefault(seqid, []).extend((length, action, start, end) for start, end in ranges)
    return adjustments


def _kept_window(seqid, trims):
    """Return the ``(first, last)`` original coordinates of ``seqid`` that remain, or None if nothing does."""
    first, last = 1, float("inf")
    for length, action, start, end in trims:
        if action == "ACTION_EXCLUDE" or (start == 1 and end == length):
            return None
        if start == 1:
            first = max(first, end + 1)
        elif end == length:
            last = min(last, start - 1)
        else:
            logger.warning(f"{seqid}:{start}..{end} ({action}) is inside the contig; features there are not adjusted")
    return (first, last) if first <= last else None


def _clip_feature(feat, first, last):
    """Return ``feat`` clipped to ``first..last`` and shifted so ``first`` becomes 1, or None if nothing is left.

    The first coordinate is the feature's 5' end on either strand, so it gets ``<``
    when that end was cut off, and the last coordinate gets ``>`` when the 3' end was.
    """
    offset = first - 1
    kept = []
    for i, (start_raw, end_raw) in enumerate(feat.intervals):
        start, end = int(start_raw.lstrip("<>")), int(end_raw.lstrip("<>"))
        low, high = sorted((start, end))
        if high < first or low > last:
            continue
        low, high = max(low, first), min(high, last)
        new_start, new_end = (high, low) if start > end else (low, high)
        kept.append((i, new_start, new_end, new_start != start, new_end != end))
    if not kept:
        return None

    intervals = [[str(start - offset), str(end - offset)] for _, start, end, _, _ in kept]
    first_idx, _, _, start_cut, _ = kept[0]
    last_idx, _, _, _, end_cut = kept[-1]
    start_marker = feat.intervals[0][0][0] if feat.intervals[0][0][0] in "<>" else ""
    end_marker = feat.intervals[-1][1][0] if feat.intervals[-1][1][0] in "<>" else ""
    if first_idx != 0 or start_cut:
        start_marker = "<"
    if last_idx != len(feat.intervals) - 1 or end_cut:
        end_marker = ">"
    intervals[0][0] = start_marker + intervals[0][0]
    intervals[-1][1] = end_marker + intervals[-1][1]
    return TblFeature(feat.key, intervals, list(feat.qualifiers))
