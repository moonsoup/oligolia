"""#86 — the Restriction and Digest tabs must agree about a pasted template.

The GUI put the two side by side: the Restriction tab's Positions column comes
from `/primers/restriction_sites`, the Digest tab's Start/End columns and its
"(template: N bp)" caption from `/primers/digest`, on the same template box. A
GenBank ORIGIN paste or a CRLF sequence used to shift one relative to the other,
because the panel pre-normalised with the same two-character strip the routers
did.

Driven the way `test_design_85.py` drives the real widget: build the panel under
the offscreen platform plugin, type into the box, click the buttons, and assert
on what the tables actually say.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from gui.panels.primers_panel import PrimersPanel  # noqa: E402
from backend.tests._puc19 import PUC19  # noqa: E402

#: pUC19's single EcoRI site, and the cut one base into it (EcoRI fst5 == 1).
ECORI_SITE = 395
ECORI_CUT = 396


def genbank_origin(seq: str) -> str:
    """The sequence as pasted out of a GenBank ORIGIN block."""
    lines = []
    for i in range(0, len(seq), 60):
        blocks = " ".join(seq[i + j:i + j + 10].lower() for j in range(0, 60, 10))
        lines.append(f"{i + 1:>9} {blocks}")
    return "\n".join(lines)


def crlf(seq: str) -> str:
    return "\r\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


@pytest.fixture
def panel(app: QApplication) -> PrimersPanel:
    return PrimersPanel()


def _cell(table, row: int, col: int) -> str:
    item = table.item(row, col)
    return item.text() if item is not None else ""


def _ecori_row(table) -> int:
    for row in range(table.rowCount()):
        if _cell(table, row, 0) == "EcoRI":
            return row
    raise AssertionError("EcoRI is not in the Restriction table")


@pytest.mark.parametrize("flavour", ["clean", "genbank_origin", "crlf"])
def test_restriction_and_digest_tabs_agree_on_a_pasted_template(panel, flavour) -> None:
    template = {"clean": PUC19,
                "genbank_origin": genbank_origin(PUC19),
                "crlf": crlf(PUC19)}[flavour]
    panel._template.setPlainText(template)
    panel._dig_enzymes.setCurrentText("EcoRI")

    panel._run_restriction()
    panel._run_digest()

    # Restriction tab: the site is where it is in the molecule.
    assert _cell(panel._re_table, _ecori_row(panel._re_table), 3) == str(ECORI_SITE)

    # Digest tab: the caption is the molecule's length, and the fragment that
    # begins at the cut begins at the cut the Restriction tab implies.
    assert f"{len(PUC19):,} bp" in panel._dig_status.text()
    starts = {_cell(panel._dig_table, r, 2) for r in range(panel._dig_table.rowCount())}
    assert str(ECORI_SITE + 1) == str(ECORI_CUT)
    assert str(ECORI_CUT) in starts
    sizes = sorted(int(_cell(panel._dig_table, r, 1).replace(",", ""))
                   for r in range(panel._dig_table.rowCount()))
    assert sizes == [396, 2290]


def test_digest_fragments_handed_to_assembly_are_dna(panel) -> None:
    """`_last_digest` is what the Assembly tab ligates — it must be sequence."""
    panel._template.setPlainText(genbank_origin(PUC19))
    panel._dig_enzymes.setCurrentText("EcoRI")
    panel._run_digest()

    assert panel._last_digest is not None
    for frag in panel._last_digest.fragments:
        assert sorted({c for c in frag.sequence if c not in "ACGT"}) == []


def test_assembly_fragment_lines_are_normalised_per_line(panel) -> None:
    """One fragment per line still, with the layout characters gone from each."""
    panel._asm_input.setPlainText("  1 acgt acgt\r\n\r\n 61 ttttGGGG\n   \n")
    assert panel._asm_lines() == ["ACGTACGT", "TTTTGGGG"]
