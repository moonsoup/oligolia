"""#95: through the panel, a point mutation inside a CDS must not be silent.

`_commit_edit` is the only place the shipped app edits a sequence. It used to
hand `shift_annotations` the coordinates and nothing else, so a CDS that
contained the edit came back with its coordinates adjusted and its pre-edit
`/translation` intact — the exported record declaring a protein its own bases no
longer encode. Handing over the edited bases lets the claim be restated instead,
and the user is told it was.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from backend.models.sequence import Annotation, MoleculeType, Sequence, Strand  # noqa: E402
from gui.panels import sequence_panel as sp  # noqa: E402

CDS_BASES = "ATG" + "GCT" * 10 + "TAA"
PROTEIN = "M" + "A" * 10
SEQ = "CCCCC" + CDS_BASES + "GGGGG"
CDS_START, CDS_END = 5, 5 + len(CDS_BASES)


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


def _seq() -> Sequence:
    return Sequence(
        id="demo", name="demo", seq=SEQ, length=len(SEQ), molecule_type=MoleculeType.DNA,
        annotations=[
            Annotation(feature_type="CDS", start=CDS_START, end=CDS_END,
                       parts=[(CDS_START, CDS_END)], strand=Strand.PLUS,
                       qualifiers={"gene": "orf", "translation": PROTEIN,
                                   "codon_start": "1"}),
        ],
    )


def _panel(app, sequence: Sequence) -> sp.SequencePanel:
    panel = sp.SequencePanel()
    panel.add_sequence(sequence)
    panel._active = sequence
    panel._render_active()
    return panel


def _choose(panel, op: str) -> None:
    panel._op_combo.setCurrentIndex(
        next(i for i in range(panel._op_combo.count())
             if panel._op_combo.itemData(i) == op)
    )


def test_a_codon_swap_restates_the_translation_and_says_so(app, tmp_path) -> None:
    seq = _seq()
    panel = _panel(app, seq)

    _choose(panel, "replace")
    panel._pos_start.setValue(CDS_START + 3)   # the second codon
    panel._pos_end.setValue(CDS_START + 6)
    panel._insert_seq.setText("TGG")
    panel._apply_op()

    assert len(seq.annotations) == 1, "the CDS was dropped instead of restated"
    cds = seq.annotations[0]
    assert (cds.start, cds.end) == (CDS_START, CDS_END), "an equal-length edit moved it"
    assert cds.qualifiers["translation"] == "MW" + "A" * 9

    shown = panel._result_display.toPlainText()
    assert "restated" in shown and "orf" in shown, shown

    shot = tmp_path / "restated.png"
    panel.grab().save(str(shot))
    assert shot.stat().st_size > 0


def test_an_edit_outside_the_cds_leaves_its_translation_alone(app) -> None:
    seq = _seq()
    panel = _panel(app, seq)

    _choose(panel, "insert")
    panel._pos_start.setValue(0)
    panel._insert_seq.setText("TTTTTT")
    panel._apply_op()

    cds = seq.annotations[0]
    assert (cds.start, cds.end) == (CDS_START + 6, CDS_END + 6)
    assert cds.qualifiers["translation"] == PROTEIN
    assert "restated" not in panel._result_display.toPlainText()
