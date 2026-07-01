"""Tests for the shared feature-type color map (gui/panels/feature_colors.py)."""

from backend.models.sequence import Annotation
from gui.panels.feature_colors import FEATURE_PALETTE, feature_color_map


def _ann(ftype: str) -> Annotation:
    return Annotation(feature_type=ftype, start=0, end=3)


def test_empty_annotations_give_empty_map() -> None:
    assert feature_color_map([]) == {}


def test_one_color_per_distinct_type() -> None:
    anns = [_ann("CDS"), _ann("CDS"), _ann("exon"), _ann("promoter")]
    m = feature_color_map(anns)
    assert set(m.keys()) == {"CDS", "exon", "promoter"}


def test_assignment_is_deterministic_regardless_of_order() -> None:
    a = feature_color_map([_ann("exon"), _ann("CDS")])
    b = feature_color_map([_ann("CDS"), _ann("exon")])
    assert {k: v.name() for k, v in a.items()} == {k: v.name() for k, v in b.items()}


def test_colors_come_from_palette_and_cycle() -> None:
    anns = [_ann(f"type{i:02d}") for i in range(len(FEATURE_PALETTE) + 2)]
    m = feature_color_map(anns)
    palette = {c.lower() for c in FEATURE_PALETTE}
    assert all(color.name().lower() in palette for color in m.values())
    # Wrap-around: the 9th sorted type reuses the first palette entry.
    types = sorted(a.feature_type for a in anns)
    assert m[types[0]].name().lower() == m[types[len(FEATURE_PALETTE)]].name().lower()
