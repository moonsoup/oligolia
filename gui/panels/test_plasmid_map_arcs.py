"""#57: a spliced or origin-spanning feature must not be drawn as one arc.

`join(91..100,1..20)` on a 100 bp circular record collapsed to outer bounds
[0,100) and was drawn as a full circle. With Annotation.parts available, the map
draws one arc per part; `arc_spans` is the pure geometry that decides it.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from backend.models.sequence import Annotation, Strand  # noqa: E402
from gui.panels.plasmid_map import PlasmidMapWidget  # noqa: E402


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def widget(app: QApplication) -> PlasmidMapWidget:
    w = PlasmidMapWidget()
    w._length = 100
    return w


def test_a_simple_feature_is_one_arc(widget: PlasmidMapWidget) -> None:
    ann = Annotation(feature_type="gene", start=10, end=30, parts=[(10, 30)], strand=Strand.PLUS)
    assert widget.arc_spans(ann) == [(10, 20)]


def test_a_feature_with_no_parts_falls_back_to_its_bounds(widget: PlasmidMapWidget) -> None:
    """Anything constructed before `parts` existed must still draw."""
    ann = Annotation(feature_type="gene", start=10, end=30, strand=Strand.PLUS)
    assert widget.arc_spans(ann) == [(10, 20)]


def test_a_spliced_feature_draws_its_exons_not_the_intron(widget: PlasmidMapWidget) -> None:
    ann = Annotation(feature_type="CDS", start=4, end=40, parts=[(4, 10), (29, 40)], strand=Strand.PLUS)
    assert widget.arc_spans(ann) == [(4, 6), (29, 11)]


def test_an_origin_spanning_feature_is_two_arcs_not_a_full_circle(widget: PlasmidMapWidget) -> None:
    """The exact #57 case."""
    ann = Annotation(feature_type="gene", start=0, end=100, parts=[(90, 100), (0, 20)], strand=Strand.PLUS)
    spans = widget.arc_spans(ann)

    assert spans == [(90, 10), (0, 20)], spans
    total = sum(span for _start, span in spans)
    assert total == 30, f"the feature covers 30 bp, drawn as {total}"
    assert total < widget._length, "must not draw the whole plasmid"


def test_a_zero_length_part_still_draws_something(widget: PlasmidMapWidget) -> None:
    """Degenerate but real in some files; must not vanish or divide by zero."""
    ann = Annotation(feature_type="misc", start=5, end=5, parts=[(5, 5)], strand=Strand.PLUS)
    assert widget.arc_spans(ann) == [(5, 1)]
