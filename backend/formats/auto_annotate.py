"""Homology-based auto-annotation of common vector parts (issue #43).

Scans a sequence against the bundled reference-parts library and returns
candidate :class:`Annotation`s for near-matches on both strands. v1 uses a
simple mismatch-tolerant (Hamming) sliding window — no BLAST dependency —
which is ample for plasmid-sized input (a few kb to ~20 kb). For genome-scale
input a BLAST/k-mer-index upgrade would matter; noted here as future work.
"""

from __future__ import annotations

from ..models.sequence import Annotation, ReferenceFeature, Strand
from .features_library import load_common_features

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _best_windows(haystack: str, needle: str, max_mm: int) -> list[tuple[int, int]]:
    """(start, mismatches) for every window of ``needle`` within ``max_mm``."""
    n, m = len(haystack), len(needle)
    hits: list[tuple[int, int]] = []
    for i in range(n - m + 1):
        mm = 0
        for a, b in zip(haystack[i:i + m], needle):
            if a != b:
                mm += 1
                if mm > max_mm:
                    break
        if mm <= max_mm:
            hits.append((i, mm))
    return hits


def auto_annotate(
    seq: str,
    library: list[ReferenceFeature] | None = None,
    *,
    min_identity: float = 0.90,
    min_length: int = 12,
) -> list[Annotation]:
    """Return candidate annotations for library parts found in ``seq``.

    Matches are searched on both strands with up to ``(1 - min_identity)`` of
    each part's length allowed as mismatches. Parts shorter than
    ``min_length`` are skipped to avoid spurious short hits. Results carry
    ``qualifiers["auto_detected"] = "true"`` so callers can distinguish them
    from user-added annotations.
    """
    if library is None:
        library = load_common_features()
    seq = seq.upper().replace(" ", "").replace("\n", "")
    n = len(seq)
    rc = _reverse_complement(seq)

    found: list[Annotation] = []
    seen: set[tuple[int, int, str, str]] = set()

    for part in library:
        needle = part.sequence.upper()
        m = len(needle)
        if m < min_length or m > n:
            continue
        max_mm = int(m * (1.0 - min_identity))

        # Forward strand.
        for start, mm in _best_windows(seq, needle, max_mm):
            key = (start, start + m, part.name, "+")
            if key in seen:
                continue
            seen.add(key)
            found.append(_make_annotation(part, start, start + m, Strand.PLUS, m, mm))

        # Reverse strand: map rc coordinates back to forward orientation.
        for rstart, mm in _best_windows(rc, needle, max_mm):
            fstart = n - (rstart + m)
            key = (fstart, fstart + m, part.name, "-")
            if key in seen:
                continue
            seen.add(key)
            found.append(_make_annotation(part, fstart, fstart + m, Strand.MINUS, m, mm))

    found.sort(key=lambda a: (a.start, a.end))
    return found


def _make_annotation(
    part: ReferenceFeature, start: int, end: int, strand: Strand, length: int, mm: int
) -> Annotation:
    identity = round((length - mm) / length, 3)
    return Annotation(
        feature_type=part.feature_type,
        start=start,
        end=end,
        strand=strand,
        qualifiers={
            "auto_detected": "true",
            "label": part.name,
            "identity": str(identity),
            "mismatches": str(mm),
            "source": part.source,
        },
    )
