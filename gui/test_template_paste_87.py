"""#87 — a pasted FASTA header must produce a readable refusal, in both tabs.

The Digest tab used to show `Error: Invalid character found in >PUC19TCGCGC…` —
a Biopython internal leaking out of an unhandled TypeError — while the
Restriction tab, on the very same box, filled its table with EcoRI at 401,
counting `>pUC19`'s letters as bases. Now both routers raise HTTPException(400)
from one shared guard, so the panel has to surface that `detail` and leave no
table of wrong numbers behind.

Driven like `test_template_paste_86.py`: the real widget under the offscreen
platform plugin, clicked, then asserted on what the status lines actually say.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from gui.panels.primers_panel import PrimersPanel  # noqa: E402
from backend.tests._puc19 import PUC19  # noqa: E402

FASTA_PASTE = ">pUC19\n" + PUC19
#: Where the Restriction tab reported EcoRI while counting the header as bases.
ECORI_SITE_OVER_HEADER = 401


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


@pytest.fixture
def panel(app: QApplication) -> PrimersPanel:
    return PrimersPanel()


def _table_text(table) -> str:
    return " ".join(
        table.item(r, c).text() if table.item(r, c) is not None else ""
        for r in range(table.rowCount())
        for c in range(table.columnCount())
    )


def test_both_tabs_explain_the_refusal_rather_than_leaking_biopython(panel) -> None:
    panel._template.setPlainText(FASTA_PASTE)
    panel._dig_enzymes.setCurrentText("EcoRI")

    panel._run_restriction()
    panel._run_digest()

    for name, status in (("Restriction", panel._re_status), ("Digest", panel._dig_status)):
        text = status.text()
        assert "Invalid character found in" not in text, (name, text)
        assert "'>'" in text, (name, text)
        assert "FASTA" in text, (name, text)
        # HTTPException's own stringification would prefix the status code.
        assert "400:" not in text, (name, text)


def test_neither_tab_shows_numbers_counted_over_the_header(panel) -> None:
    """A plausible-looking table of wrong offsets is worse than an error."""
    panel._template.setPlainText(FASTA_PASTE)
    panel._dig_enzymes.setCurrentText("EcoRI")

    panel._run_restriction()
    panel._run_digest()

    assert panel._re_table.rowCount() == 0
    assert panel._dig_table.rowCount() == 0
    assert str(ECORI_SITE_OVER_HEADER) not in _table_text(panel._re_table)


def test_a_stale_good_result_is_cleared_by_the_refusal(panel) -> None:
    """The user edits the box and re-runs; the old rows must not survive."""
    panel._template.setPlainText(PUC19)
    panel._dig_enzymes.setCurrentText("EcoRI")
    panel._run_restriction()
    panel._run_digest()
    assert panel._re_table.rowCount() > 0
    assert panel._dig_table.rowCount() > 0

    panel._template.setPlainText(FASTA_PASTE)
    panel._run_restriction()
    panel._run_digest()
    assert panel._re_table.rowCount() == 0
    assert panel._dig_table.rowCount() == 0


def test_a_layout_heavy_paste_still_works_in_both_tabs(panel) -> None:
    """The #86 path is untouched: CRLF wrapping is normalised, not refused."""
    panel._template.setPlainText("\r\n".join(PUC19[i:i + 60] for i in range(0, len(PUC19), 60)))
    panel._dig_enzymes.setCurrentText("EcoRI")

    panel._run_restriction()
    panel._run_digest()

    assert "Error" not in panel._re_status.text()
    assert "Error" not in panel._dig_status.text()
    assert panel._re_table.rowCount() > 0
    assert f"{len(PUC19):,} bp" in panel._dig_status.text()
