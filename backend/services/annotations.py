"""Keep feature coordinates correct when a sequence is edited.

#59: insert, delete, replace and reverse-complement all left annotations exactly
where they were, so after any edit every feature pointed at the wrong bases — and
export wrote them that way. A GenBank file round-tripped through one insertion
came out describing a different molecule.

Pure and free of Qt, so the arithmetic is testable on its own and any other edit
path can reuse it.

Every function works on `parts` (from #57), not just the outer bounds, because a
spliced feature's exons move independently: an edit inside the intron shifts the
downstream exon and leaves the upstream one alone.
"""

from __future__ import annotations

from ..models.sequence import Annotation, Strand


def _parts_of(ann: Annotation) -> list[tuple[int, int]]:
    return list(ann.parts) if ann.parts else [(ann.start, ann.end)]


def _rebuilt(ann: Annotation, parts: list[tuple[int, int]], strand: Strand | None = None) -> Annotation:
    """A copy of `ann` with new parts, and outer bounds derived from them."""
    return ann.model_copy(update={
        "parts": parts,
        "start": min(s for s, _e in parts),
        "end": max(e for _s, e in parts),
        "strand": strand if strand is not None else ann.strand,
    })


def shift_annotations(
    annotations: list[Annotation],
    *,
    start: int,
    end: int,
    inserted: int,
) -> tuple[list[Annotation], list[Annotation]]:
    """Move annotations across an edit that replaces `[start, end)` with `inserted` bases.

    Insert is `start == end`; delete is `inserted == 0`; replace is the general
    case. Returns `(kept, dropped)`.

    A feature whose bases were themselves edited is **dropped**, not truncated.
    Truncating would keep a feature whose sequence no longer matches what it
    claims to describe, which is the same class of quiet wrongness as the
    original bug — so the caller gets the list and can say which features went.

    A feature that merely *contains* the edit grows or shrinks with it: those
    bases are still its bases. A feature strictly after the edit shifts by the
    net length change; one strictly before it does not move.
    """
    if end < start:
        raise ValueError(f"edit end ({end}) is before start ({start})")

    delta = inserted - (end - start)
    kept: list[Annotation] = []
    dropped: list[Annotation] = []

    for ann in annotations:
        new_parts: list[tuple[int, int]] = []
        lost = False

        for p_start, p_end in _parts_of(ann):
            if p_end <= start:
                # Entirely before the edit.
                new_parts.append((p_start, p_end))
            elif p_start >= end:
                # Entirely after it.
                new_parts.append((p_start + delta, p_end + delta))
            elif p_start <= start and p_end >= end:
                # Contains the edit: keep the start, move the end.
                new_parts.append((p_start, p_end + delta))
            else:
                # Overlaps the edited region: its own bases changed.
                lost = True
                break

        if lost:
            dropped.append(ann)
        else:
            kept.append(_rebuilt(ann, new_parts))

    return kept, dropped


def flip_annotations(annotations: list[Annotation], seq_len: int) -> list[Annotation]:
    """Map annotations onto the reverse complement of a `seq_len`-long sequence.

    A half-open interval `[s, e)` on the forward strand becomes
    `[seq_len - e, seq_len - s)` on the reverse, the strand inverts, and a spliced
    feature's parts come back in the opposite order — what was the first exon is
    the last one in the new coordinates.

    Applying this twice returns the original, which is how the mapping is checked
    rather than reasoned about.
    """
    flipped: list[Annotation] = []

    for ann in annotations:
        parts = [(seq_len - p_end, seq_len - p_start) for p_start, p_end in _parts_of(ann)]
        parts.reverse()
        if ann.strand == Strand.PLUS:
            strand = Strand.MINUS
        elif ann.strand == Strand.MINUS:
            strand = Strand.PLUS
        else:
            strand = ann.strand  # BOTH / unstranded has no orientation to flip
        flipped.append(_rebuilt(ann, parts, strand))

    return flipped
