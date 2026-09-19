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

from typing import NamedTuple

from Bio.Seq import Seq

from ..models.sequence import Annotation, LocationPart, Strand

#: `<` on a lower boundary becomes `>` on the upper one when the strand is
#: reflected, and vice versa — the same mapping Biopython's `BeforePosition._flip`
#: makes (#94).
_FLIPPED_CLASS = {"before": "after", "after": "before", "exact": "exact"}


def _parts_of(ann: Annotation) -> list[tuple[int, int]]:
    return list(ann.parts) if ann.parts else [(ann.start, ann.end)]


def _rebuilt(
    ann: Annotation,
    parts: list[tuple[int, int]],
    strand: Strand | None = None,
    part_details: list[LocationPart] | None = None,
) -> Annotation:
    """A copy of `ann` with new parts, and outer bounds derived from them."""
    update: dict = {
        "parts": parts,
        "start": min(s for s, _e in parts),
        "end": max(e for _s, e in parts),
        "strand": strand if strand is not None else ann.strand,
    }
    if part_details is not None:
        update["part_details"] = part_details
    return ann.model_copy(update=update)


def extracted_bases(ann: Annotation, sequence: str) -> str:
    """The subsequence `ann` denotes in `sequence`.

    Parts are concatenated in the order they are stored — the feature's own
    5'-to-3' reading order (#93) — and each is reverse-complemented on the minus
    strand, which is what `Bio.SeqFeature.extract` does and therefore what the
    exported record will mean.
    """
    pieces = [sequence[p_start:p_end] for p_start, p_end in _parts_of(ann)]
    if ann.strand == Strand.MINUS:
        return "".join(str(Seq(piece).reverse_complement()) for piece in pieces)
    return "".join(pieces)


def restated_translation(ann: Annotation, sequence: str) -> str | None:
    """The protein `ann`'s bases in `sequence` actually encode, or None.

    None means the claim cannot be restated — there is nothing left to
    translate — and the caller should drop the feature rather than carry a
    protein that no longer follows from the bases under it (#95).

    `/codon_start` and `/transl_table` are read from the feature itself, so a
    record that declared a shifted frame or a non-standard code keeps it. A
    trailing partial codon is dropped before translating, which is what
    `Bio.Seq.translate` does with it anyway, minus the warning.
    """
    bases = extracted_bases(ann, sequence)
    try:
        codon_start = int(str(ann.qualifiers.get("codon_start", 1)))
    except (TypeError, ValueError):
        codon_start = 1
    try:
        table = int(str(ann.qualifiers.get("transl_table", 1)))
    except (TypeError, ValueError):
        table = 1

    coding = bases[max(codon_start - 1, 0):]
    coding = coding[: len(coding) // 3 * 3]
    if not coding:
        return None
    return str(Seq(coding).translate(table=table, to_stop=True))


def shift_annotations(
    annotations: list[Annotation],
    *,
    start: int,
    end: int,
    inserted: int,
    new_sequence: str | None = None,
) -> tuple[list[Annotation], list[Annotation]]:
    """Move annotations across an edit that replaces `[start, end)` with `inserted` bases.

    Insert is `start == end`; delete is `inserted == 0`; replace is the general
    case. `new_sequence` is the sequence *after* the edit; see below for the one
    thing it changes. Returns `(kept, dropped)`.

    A feature whose bases were themselves edited is **dropped**, not truncated.
    Truncating would keep a feature whose sequence no longer matches what it
    claims to describe, which is the same class of quiet wrongness as the
    original bug — so the caller gets the list and can say which features went.

    A feature that merely *contains* the edit grows or shrinks with it: those
    bases are still its bases. A feature strictly after the edit shifts by the
    net length change; one strictly before it does not move.

    Two ways that containing branch used to keep a feature it should not have
    (#95):

    * **Nothing left to describe.** Deleting exactly a feature's own bases
      satisfies `p_start <= start and p_end >= end`, so the feature never reached
      the overlap branch: it was kept, collapsed to zero length, and exported as
      the INSDC between-position `5073^5074` — a site *between* two bases, which
      the record never claimed — or, at the very start of a record, as the
      impossible `0^1`. A part the edit empties is gone, and a feature with no
      part left is dropped like any other feature whose bases were edited.
    * **A claim about the exact bases.** Coordinates are only part of what a
      feature asserts. `/translation` states the protein *these* bases encode, so
      an edit inside a CDS falsifies it however the interval is adjusted — an
      equal-length replacement inside a CDS moves no coordinate at all and still
      leaves the record declaring the pre-edit protein. Given `new_sequence` the
      translation is restated from the edited bases; without it there is no way
      to restate it, so the feature is dropped and reported rather than kept
      saying something false.

    Features that merely surround the edit without any of their own bases
    changing — an insert flush against a boundary — are untouched by both rules,
    as are features that make no base-level claim.
    """
    if end < start:
        raise ValueError(f"edit end ({end}) is before start ({start})")

    delta = inserted - (end - start)
    kept: list[Annotation] = []
    dropped: list[Annotation] = []

    for ann in annotations:
        new_parts: list[tuple[int, int]] = []
        lost = False
        bases_changed = False

        for p_start, p_end in _parts_of(ann):
            if p_end <= start:
                # Entirely before the edit.
                new_parts.append((p_start, p_end))
            elif p_start >= end:
                # Entirely after it.
                new_parts.append((p_start + delta, p_end + delta))
            elif p_start <= start and p_end >= end:
                # Contains the edit: keep the start, move the end. Getting here
                # means the edit is strictly interior — the two branches above
                # already took everything flush with a boundary — so these bases
                # did change, whatever the arithmetic does to the interval.
                shifted = (p_start, p_end + delta)
                if end > start and shifted[1] <= shifted[0]:
                    # The edit removed every base this part described.
                    lost = True
                    break
                new_parts.append(shifted)
                bases_changed = True
            else:
                # Overlaps the edited region: its own bases changed.
                lost = True
                break

        if lost:
            dropped.append(ann)
            continue

        rebuilt = _rebuilt(ann, new_parts)
        if bases_changed and ann.qualifiers.get("translation"):
            protein = restated_translation(rebuilt, new_sequence) if new_sequence else None
            if protein is None:
                dropped.append(ann)
                continue
            rebuilt = rebuilt.model_copy(
                update={"qualifiers": {**ann.qualifiers, "translation": protein}}
            )
        kept.append(rebuilt)

    return kept, dropped


class SplicedAnnotations(NamedTuple):
    """What a splice did to a record's annotations.

    `kept` and `dropped` are `shift_annotations`' own two lists. `restated` holds
    the *original* objects whose `/translation` was rewritten from the edited
    bases (#95), so a caller can name them; it is a subset of the features
    `kept` corresponds to, in the same order.
    """

    kept: list[Annotation]
    dropped: list[Annotation]
    restated: list[Annotation]


def spliced_annotations(
    annotations: list[Annotation],
    *,
    start: int,
    end: int,
    inserted: int,
    new_sequence: str | None = None,
) -> SplicedAnnotations:
    """`shift_annotations`, plus which features had their translation restated.

    The GUI's `_commit_edit` worked this out inline, so the HTTP edit endpoint
    had nothing to call and built its stored record from five fields instead —
    dropping every annotation and the molecule's topology on the way (#96).
    Both paths call this now, so an edit means the same thing whichever one of
    them the user reached it through.
    """
    kept, dropped = shift_annotations(
        annotations, start=start, end=end, inserted=inserted, new_sequence=new_sequence
    )
    # `dropped` holds the original objects and both lists keep their input
    # order, so the survivors line up one-for-one with `kept`.
    lost_ids = {id(a) for a in dropped}
    survivors = [a for a in annotations if id(a) not in lost_ids]
    restated = [
        a for a, k in zip(survivors, kept)
        if a.qualifiers.get("translation") != k.qualifiers.get("translation")
    ]
    return SplicedAnnotations(kept, dropped, restated)


def flip_annotations(annotations: list[Annotation], seq_len: int) -> list[Annotation]:
    """Map annotations onto the reverse complement of a `seq_len`-long sequence.

    A half-open interval `[s, e)` on the forward strand becomes
    `[seq_len - e, seq_len - s)` on the reverse, and the strand inverts.

    The **part list keeps its order** (#93). A `CompoundLocation`'s parts are
    stored in the feature's own 5'-to-3' reading order, not in ascending
    coordinate order — Biopython's INSDC writer spells the convention out: for a
    minus-strand join it "expect[s] the CompoundLocation and its parts to all be
    marked as strand == -1, and to be in the order 19:100 then 0:10". Reflecting
    each interval through `seq_len - x` already preserves that reading order, so
    reversing the list afterwards undid it: a spliced CDS came back with its exons
    swapped, which is the same length, the same bases and the same strand, but a
    different protein — and it reached the saved file with nothing to signal it.

    Applying this twice returns the original, but that round trip is *not* what
    checks the mapping: reversing a list twice restores it whether or not
    reversing it was right in the first place. What pins it is that every feature
    still extracts the same string, which is what
    `backend/tests/test_reverse_complement_part_order_93.py` asserts against
    `Bio.SeqRecord.reverse_complement(features=True)`.

    A partial boundary travels with the boundary it describes (#94): reflecting
    `[s, e)` makes the old upper boundary the new lower one, so a `>` on the end
    comes back as a `<` on the start — exactly what `SimpleLocation._flip` does.
    """
    flipped: list[Annotation] = []

    for ann in annotations:
        parts = [(seq_len - p_end, seq_len - p_start) for p_start, p_end in _parts_of(ann)]
        details = [
            LocationPart(
                start_class=_FLIPPED_CLASS.get(d.end_class, "exact"),
                end_class=_FLIPPED_CLASS.get(d.start_class, "exact"),
                ref=d.ref,
            )
            for d in ann.details_per_part()
        ]
        if ann.strand == Strand.PLUS:
            strand = Strand.MINUS
        elif ann.strand == Strand.MINUS:
            strand = Strand.PLUS
        else:
            strand = ann.strand  # BOTH / unstranded has no orientation to flip
        flipped.append(_rebuilt(ann, parts, strand, details))

    return flipped
