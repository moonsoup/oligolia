"""The primers panel says what its Tm column assumes (#81, item 5).

Driven through the real widget offscreen rather than asserted on internals: build
it in a QApplication, paste a template, click Design Primers, wait for the worker,
and read the caption that ends up next to the table.

The complaint this closes (#79.4) is that the table showed a Tm column and said
nothing about the buffer behind it — and the same primer reads about 10 degC apart
between the pinned bench buffer and a PCR-like one.
"""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from backend.routers.primers import TM_CONDITIONS  # noqa: E402
from gui.panels.primers_panel import PrimersPanel  # noqa: E402

TEMPLATE = (
    "ATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACGTGGCATCACGTTT"
    "CGGCCATGGACGTAGGCATCACGTGGCATCACGTGGCCGGCGGCGGCGCTGCTGCAAGCTT"
    "ACGTGGCATCACGATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACG"
)


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def panel(app):
    return PrimersPanel()


def test_the_tm_column_names_its_method_and_buffer_before_any_run(panel) -> None:
    """Stated up front, not only after a successful design."""
    note = panel._pcr_tm_note.text()
    assert "nearest-neighbour" in note
    assert "DNA_NN3" in note
    assert f"Na⁺ {TM_CONDITIONS['Na']:g} mM" in note
    assert f"Mg²⁺ {TM_CONDITIONS['Mg']:g} mM" in note
    assert f"primer {TM_CONDITIONS['dnac1']:g} nM" in note


def test_the_caption_reports_the_buffer_the_shown_rows_were_computed_under(panel, tmp_path) -> None:
    """End to end through the widget: design, then read the caption by the table."""
    panel._template.setPlainText(TEMPLATE)
    panel._prod_min.setValue(100)
    panel._prod_max.setValue(200)
    panel._tm_min.setValue(45.0)
    panel._tm_max.setValue(75.0)
    panel._run_pcr()

    # The design runs on a QThread; wait for it rather than sleeping.
    assert panel._worker.wait(60_000), "primer design did not finish"
    QApplication.processEvents()

    assert panel._pcr_table.rowCount() > 0, "no primer pairs reached the table"
    note = panel._pcr_tm_note.text()
    assert "nearest-neighbour" in note, note
    assert f"Mg²⁺ {TM_CONDITIONS['Mg']:g} mM" in note, note

    panel.grab().save(str(tmp_path / "primers_tm_note.png"))
