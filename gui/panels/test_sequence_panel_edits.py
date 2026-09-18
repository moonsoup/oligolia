"""#59: drive the real SequencePanel through an edit and check the features moved.

The pure arithmetic is covered in backend/tests/test_annotation_edits.py. This
checks the wiring: that the panel actually calls it, that the undo checkpoint
still works, that a bad range is refused, and that duplicate IDs no longer
collapse.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from backend.models.sequence import Annotation, MoleculeType, Sequence, Strand  # noqa: E402
from gui.panels import sequence_panel as sp  # noqa: E402


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def recorded(monkeypatch) -> list[tuple[str, str]]:
    seen: list[tuple[str, str]] = []

    class FakeBox:
        @staticmethod
        def warning(_p, title, text, *a, **k):
            seen.append((title, text))

        @staticmethod
        def critical(_p, title, text, *a, **k):
            seen.append((title, text))

        @staticmethod
        def information(_p, title, text, *a, **k):
            seen.append((title, text))

    monkeypatch.setattr(sp, "QMessageBox", FakeBox)
    return seen


SEQ = "ATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACGTGG"  # 51 nt


def _seq(seq_id: str = "demo", molecule=MoleculeType.DNA) -> Sequence:
    return Sequence(
        id=seq_id, name=seq_id, seq=SEQ, length=len(SEQ), molecule_type=molecule,
        annotations=[
            Annotation(feature_type="gene", start=30, end=40, parts=[(30, 40)],
                       strand=Strand.PLUS, qualifiers={"gene": "downstream"}),
            Annotation(feature_type="gene", start=0, end=5, parts=[(0, 5)],
                       strand=Strand.PLUS, qualifiers={"gene": "upstream"}),
        ],
    )


def _panel(app, sequence: Sequence) -> sp.SequencePanel:
    panel = sp.SequencePanel()
    panel.add_sequence(sequence)
    panel._active = sequence
    panel._render_active()
    return panel


def _by_name(sequence: Sequence, name: str):
    return next((a for a in sequence.annotations if a.qualifiers.get("gene") == name), None)


def test_an_insertion_moves_the_downstream_feature(app, recorded) -> None:
    seq = _seq()
    panel = _panel(app, seq)

    panel._op_combo.setCurrentIndex(
        next(i for i in range(panel._op_combo.count())
             if panel._op_combo.itemData(i) == "insert"))
    panel._pos_start.setValue(10)
    panel._insert_seq.setText("AAAAA")
    panel._apply_op()

    assert len(seq.seq) == len(SEQ) + 5
    assert (_by_name(seq, "downstream").start, _by_name(seq, "downstream").end) == (35, 45)
    assert (_by_name(seq, "upstream").start, _by_name(seq, "upstream").end) == (0, 5)


def test_a_deletion_over_a_feature_drops_it_and_says_so(app, recorded) -> None:
    seq = _seq()
    panel = _panel(app, seq)

    panel._op_combo.setCurrentIndex(
        next(i for i in range(panel._op_combo.count())
             if panel._op_combo.itemData(i) == "delete"))
    panel._pos_start.setValue(28)
    panel._pos_end.setValue(45)
    panel._apply_op()

    assert _by_name(seq, "downstream") is None, "the deleted feature is still there"
    assert _by_name(seq, "upstream") is not None
    shown = panel._result_display.toPlainText()
    assert "removed" in shown and "downstream" in shown, shown


def test_replace_with_end_zero_is_refused(app, recorded) -> None:
    """#59: this GREW the sequence — 20 nt became 32."""
    seq = _seq()
    panel = _panel(app, seq)

    panel._op_combo.setCurrentIndex(
        next(i for i in range(panel._op_combo.count())
             if panel._op_combo.itemData(i) == "replace"))
    panel._pos_start.setValue(10)
    panel._pos_end.setValue(0)
    panel._insert_seq.setText("GGGGGGGGGGGG")
    panel._apply_op()

    assert len(seq.seq) == len(SEQ), "the sequence changed length on a refused edit"
    assert recorded, "the user was not told why"
    assert "range" in (recorded[0][0] + recorded[0][1]).lower()


def test_reverse_complement_moves_and_flips_features(app, recorded) -> None:
    seq = _seq()
    panel = _panel(app, seq)

    panel._op_combo.setCurrentIndex(
        next(i for i in range(panel._op_combo.count())
             if panel._op_combo.itemData(i) == "reverse_complement"))
    panel._apply_op()

    up = _by_name(seq, "upstream")
    assert (up.start, up.end) == (len(SEQ) - 5, len(SEQ))
    assert up.strand == Strand.MINUS


def test_reverse_complement_is_refused_on_a_protein(app, recorded) -> None:
    seq = _seq(molecule=MoleculeType.PROTEIN)
    panel = _panel(app, seq)

    panel._op_combo.setCurrentIndex(
        next(i for i in range(panel._op_combo.count())
             if panel._op_combo.itemData(i) == "reverse_complement"))
    panel._apply_op()

    assert seq.seq == SEQ, "a nucleotide op ran on a protein sequence"
    assert recorded, "the user was not warned"
    assert "protein" in (recorded[0][0] + recorded[0][1]).lower()


def test_two_records_with_the_same_id_both_survive(app) -> None:
    """#59: the second collapsed onto the first and inherited its undo history."""
    panel = sp.SequencePanel()
    panel.add_sequence(_seq("same"))
    panel.add_sequence(_seq("same"))

    assert len(panel._sequences) == 2, panel._sequences.keys()
    assert panel._list.count() == 2


def test_an_edit_is_still_undoable(app, recorded) -> None:
    seq = _seq()
    panel = _panel(app, seq)

    panel._op_combo.setCurrentIndex(
        next(i for i in range(panel._op_combo.count())
             if panel._op_combo.itemData(i) == "insert"))
    panel._pos_start.setValue(10)
    panel._insert_seq.setText("AAAAA")
    panel._apply_op()
    assert len(seq.seq) == len(SEQ) + 5

    panel._undo()
    assert panel._active.seq == SEQ
