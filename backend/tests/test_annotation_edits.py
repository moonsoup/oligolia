"""#59: editing a sequence left every feature pointing at the wrong bases.

Insert, delete, replace and reverse-complement all left coordinates untouched,
and export then wrote them that way — a silently corrupted record.

The logic is pure and lives in the backend so it can be tested without Qt, and
reused by anything else that edits a sequence.
"""

from __future__ import annotations

import pytest

from backend.models.sequence import Annotation, Strand
from backend.services.annotations import flip_annotations, shift_annotations


def _ann(start: int, end: int, name: str = "f", parts=None, strand=Strand.PLUS) -> Annotation:
    return Annotation(
        feature_type="gene", start=start, end=end,
        parts=parts or [(start, end)], strand=strand,
        qualifiers={"gene": name},
    )


# ── insertion ────────────────────────────────────────────────────────────────

def test_a_feature_after_an_insertion_moves_by_its_length() -> None:
    kept, dropped = shift_annotations([_ann(50, 60)], start=10, end=10, inserted=5)
    assert dropped == []
    assert (kept[0].start, kept[0].end) == (55, 65)
    assert kept[0].parts == [(55, 65)]


def test_a_feature_before_an_insertion_does_not_move() -> None:
    kept, _ = shift_annotations([_ann(0, 5)], start=10, end=10, inserted=5)
    assert (kept[0].start, kept[0].end) == (0, 5)


def test_a_feature_spanning_an_insertion_grows() -> None:
    """Inserting inside a feature extends it; the bases are still its bases."""
    kept, dropped = shift_annotations([_ann(0, 20)], start=10, end=10, inserted=5)
    assert dropped == []
    assert (kept[0].start, kept[0].end) == (0, 25)


def test_an_insertion_at_a_feature_boundary_does_not_silently_extend_it() -> None:
    """An insert exactly at the end is after the feature, not inside it."""
    kept, _ = shift_annotations([_ann(0, 10)], start=10, end=10, inserted=5)
    assert (kept[0].start, kept[0].end) == (0, 10)


# ── deletion ─────────────────────────────────────────────────────────────────

def test_a_feature_after_a_deletion_moves_back() -> None:
    kept, dropped = shift_annotations([_ann(50, 60)], start=10, end=20, inserted=0)
    assert dropped == []
    assert (kept[0].start, kept[0].end) == (40, 50)


def test_a_feature_inside_a_deletion_is_dropped_and_reported() -> None:
    """Its bases are gone. Keeping it pointing anywhere would be a lie."""
    kept, dropped = shift_annotations([_ann(12, 18, "doomed")], start=10, end=20, inserted=0)
    assert kept == []
    assert len(dropped) == 1
    assert dropped[0].qualifiers["gene"] == "doomed"


def test_a_feature_partly_deleted_is_dropped_not_truncated() -> None:
    """Truncating would keep a feature whose sequence no longer matches it.

    Dropping and saying so is the honest option; a half-feature silently
    retained is the #59 defect in a new form.
    """
    kept, dropped = shift_annotations([_ann(15, 30, "clipped")], start=10, end=20, inserted=0)
    assert kept == []
    assert dropped[0].qualifiers["gene"] == "clipped"


def test_a_deletion_after_a_feature_leaves_it_alone() -> None:
    kept, _ = shift_annotations([_ann(0, 10)], start=20, end=30, inserted=0)
    assert (kept[0].start, kept[0].end) == (0, 10)


# ── replacement ──────────────────────────────────────────────────────────────

def test_a_same_length_replacement_moves_nothing() -> None:
    kept, dropped = shift_annotations([_ann(50, 60)], start=10, end=20, inserted=10)
    assert dropped == []
    assert (kept[0].start, kept[0].end) == (50, 60)


def test_a_longer_replacement_pushes_later_features_out() -> None:
    kept, _ = shift_annotations([_ann(50, 60)], start=10, end=20, inserted=30)
    assert (kept[0].start, kept[0].end) == (70, 80)


def test_a_shorter_replacement_pulls_them_back() -> None:
    kept, _ = shift_annotations([_ann(50, 60)], start=10, end=20, inserted=2)
    assert (kept[0].start, kept[0].end) == (42, 52)


# ── spliced features: every part moves independently ────────────────────────

def test_each_part_of_a_spliced_feature_shifts() -> None:
    """#57 gave us `parts`; this is what it was needed for."""
    spliced = _ann(4, 40, "cds", parts=[(4, 10), (29, 40)])
    kept, dropped = shift_annotations([spliced], start=0, end=0, inserted=100)
    assert dropped == []
    assert kept[0].parts == [(104, 110), (129, 140)]
    assert (kept[0].start, kept[0].end) == (104, 140)


def test_a_spliced_feature_losing_one_part_is_dropped() -> None:
    spliced = _ann(4, 40, "cds", parts=[(4, 10), (29, 40)])
    kept, dropped = shift_annotations([spliced], start=5, end=35, inserted=0)
    assert kept == []
    assert dropped[0].qualifiers["gene"] == "cds"


def test_a_spliced_feature_with_an_edit_between_its_parts_keeps_both() -> None:
    """The intron is not the feature; editing it shifts the later exon only."""
    spliced = _ann(4, 40, "cds", parts=[(4, 10), (29, 40)])
    kept, dropped = shift_annotations([spliced], start=15, end=15, inserted=6)
    assert dropped == []
    assert kept[0].parts == [(4, 10), (35, 46)]


# ── reverse complement ───────────────────────────────────────────────────────

def test_reverse_complement_mirrors_coordinates() -> None:
    kept = flip_annotations([_ann(10, 20)], seq_len=100)
    assert (kept[0].start, kept[0].end) == (80, 90)
    assert kept[0].parts == [(80, 90)]


def test_reverse_complement_flips_the_strand() -> None:
    plus = flip_annotations([_ann(10, 20, strand=Strand.PLUS)], seq_len=100)
    assert plus[0].strand == Strand.MINUS
    minus = flip_annotations([_ann(10, 20, strand=Strand.MINUS)], seq_len=100)
    assert minus[0].strand == Strand.PLUS


def test_reverse_complement_reverses_the_part_order() -> None:
    """After flipping, the part list runs *downwards* in coordinates (#93).

    This expectation used to read `[(60, 71), (90, 96)]` — ascending, on what is
    now a minus-strand feature — and so locked in the defect #93 reports. Parts
    are stored in the feature's own 5'-to-3' reading order, and for a minus-strand
    `CompoundLocation` that means descending coordinates: the exon this feature
    read first still reads first, it has just moved to `[90, 96)`. Reflecting each
    interval preserves that order on its own; nothing extra reverses the list.
    """
    spliced = _ann(4, 40, "cds", parts=[(4, 10), (29, 40)])
    kept = flip_annotations([spliced], seq_len=100)
    assert kept[0].parts == [(90, 96), (60, 71)]
    assert (kept[0].start, kept[0].end) == (60, 96)


def test_flipping_twice_returns_the_original() -> None:
    """The invariant that makes the mapping checkable without arithmetic by hand."""
    original = [_ann(10, 20), _ann(4, 40, "cds", parts=[(4, 10), (29, 40)])]
    twice = flip_annotations(flip_annotations(original, 100), 100)
    for before, after in zip(original, twice):
        assert (after.start, after.end) == (before.start, before.end)
        assert after.parts == before.parts
        assert after.strand == before.strand


def test_a_feature_covering_the_whole_sequence_survives_a_flip() -> None:
    kept = flip_annotations([_ann(0, 100)], seq_len=100)
    assert (kept[0].start, kept[0].end) == (0, 100)


# ── guards ───────────────────────────────────────────────────────────────────

def test_an_empty_annotation_list_is_fine() -> None:
    assert shift_annotations([], start=0, end=5, inserted=0) == ([], [])
    assert flip_annotations([], seq_len=10) == []


def test_a_zero_length_edit_changes_nothing() -> None:
    kept, dropped = shift_annotations([_ann(10, 20)], start=5, end=5, inserted=0)
    assert dropped == []
    assert (kept[0].start, kept[0].end) == (10, 20)


def test_an_end_before_start_is_refused() -> None:
    """`Replace` with End=0 grew the sequence; the coordinates were never checked."""
    with pytest.raises(ValueError):
        shift_annotations([_ann(10, 20)], start=10, end=0, inserted=5)
