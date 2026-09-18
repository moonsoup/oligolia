"""#66.2: a stale UniProt accession made the panel show another protein.

`set_target` only wrote the accession field when the new value was non-empty, and
`_run_predict` reads that field — so loading a sequence with no accession of its
own looked it up under the PREVIOUS protein's accession and displayed that
protein's structure.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from gui.panels.structure_panel import StructurePanel  # noqa: E402

P53 = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAA"
OTHER = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVL"


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


def test_a_target_with_no_accession_clears_the_previous_one(app) -> None:
    panel = StructurePanel()
    panel.set_target(P53, gene_symbol="TP53", uniprot_id="P04637")
    assert panel._uniprot_input.text() == "P04637"

    panel.set_target(OTHER)

    assert panel._uniprot_input.text() == "", (
        "the previous protein's accession is still in the box; _run_predict reads "
        "it, so this sequence would be looked up as P04637 (#66)"
    )
    assert panel._gene_input.text() == ""


def test_a_new_accession_replaces_the_old_one(app) -> None:
    panel = StructurePanel()
    panel.set_target(P53, uniprot_id="P04637")
    panel.set_target(OTHER, uniprot_id="P68871")
    assert panel._uniprot_input.text() == "P68871"


def test_changing_target_clears_the_previous_result(app) -> None:
    """The structure and interaction points belong to the old protein."""
    panel = StructurePanel()
    panel.set_target(P53, uniprot_id="P04637")
    panel._result = object()          # stand in for a fetched StructureResult
    panel._points = [{"residue_index": 1}]
    panel._table.setRowCount(3)
    panel._btn_save.setEnabled(True)

    panel.set_target(OTHER)

    assert panel._result is None
    assert panel._points == []
    assert panel._table.rowCount() == 0
    assert not panel._btn_save.isEnabled()


def test_the_request_built_for_a_bare_sequence_carries_no_accession(app) -> None:
    """End of the chain: what _run_predict would actually send."""
    panel = StructurePanel()
    panel.set_target(P53, gene_symbol="TP53", uniprot_id="P04637")
    panel.set_target(OTHER)

    assert (panel._uniprot_input.text().strip() or None) is None
    assert (panel._gene_input.text().strip() or None) is None
