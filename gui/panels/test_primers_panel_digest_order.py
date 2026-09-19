"""The Digest tab answers the same whichever order the enzymes are typed (#88).

Driven through the real widget offscreen rather than asserted on internals: build
it in a QApplication, paste a template, type the enzyme box the way a user does
(the box splits on both ',' and '+'), click Digest, and read the table and the
status caption that come out.

The complaint this closes is that 'AvaI+KpnI' and 'KpnI+AvaI' — both natural
things to type — reported different end chemistry for the identical cut, and the
panel said nothing about it either way.
"""

from __future__ import annotations

import os
import sys

import pytest
from PyQt6.QtWidgets import QApplication

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from gui.panels.primers_panel import PrimersPanel  # noqa: E402

# KpnI (G_GTAC^C) and AvaI (C^YCGRG) both cut at 15 — pUC19's MCS arrangement.
SHARED_CUT = "AAAAAAAAAA" "GGTACCCGGG" "TTTTTTTTTT"


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def panel(app):
    return PrimersPanel()


def _digest(panel, enzyme_text: str) -> tuple[list[list[str]], str]:
    """Run the Digest tab as typed; return (table rows, status text)."""
    panel._template.setPlainText(SHARED_CUT)
    panel._dig_enzymes.setCurrentText(enzyme_text)
    panel._run_digest()
    QApplication.processEvents()
    rows = [
        [panel._dig_table.item(r, c).text() for c in range(panel._dig_table.columnCount())]
        for r in range(panel._dig_table.rowCount())
    ]
    return rows, panel._dig_status.text()


def test_typing_the_enzymes_the_other_way_round_gives_the_same_digest(panel, tmp_path):
    forward_rows, forward_status = _digest(panel, "AvaI+KpnI")
    forward_frags = [f.model_dump() for f in panel._last_digest.fragments]

    reverse_rows, reverse_status = _digest(panel, "KpnI+AvaI")
    reverse_frags = [f.model_dump() for f in panel._last_digest.fragments]

    assert forward_rows == reverse_rows
    assert forward_status == reverse_status
    # Including the end chemistry the table does not show but the Assembly tab
    # reads straight off _last_digest.
    assert forward_frags == reverse_frags

    panel.grab().save(str(tmp_path / "digest_enzyme_order.png"))


def test_the_panel_says_when_two_enzymes_share_a_cut(panel):
    _rows, status = _digest(panel, "AvaI+KpnI")
    assert "ambiguous end at 15" in status, status

    # And stays quiet when there is nothing ambiguous to report.
    _rows, plain = _digest(panel, "KpnI")
    assert "ambiguous" not in plain, plain
