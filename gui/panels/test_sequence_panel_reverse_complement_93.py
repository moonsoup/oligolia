"""#93: Reverse Complement through the shipped widget must not change the protein.

`backend/tests/test_reverse_complement_part_order_93.py` pins the arithmetic. This
drives the one place in the app that reaches it — `Sequence` tab → Edit Operation
→ Apply (`_apply_op` → `_commit_edit(flip=True)`) — and then saves the record the
way "Save all sequences as GenBank…" does (`gui/main_window.py` → `write_genbank`),
so the assertion is on the file a user would actually end up with.

The companion check that the feature table still lists every row passes with the
defect present, which is the point: the app gave the user no signal at all.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from io import StringIO  # noqa: E402

import pytest  # noqa: E402
from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from backend.formats.genbank import read_genbank, write_genbank  # noqa: E402
from backend.tests._spliced_circular_93 import (  # noqa: E402
    ORIGIN_SPANNING,
    SPLICED_CDS,
    SPLICED_CIRCULAR,
    SPLICED_PROTEIN,
)
from gui.panels.sequence_panel import SequencePanel  # noqa: E402


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def panel(app):
    p = SequencePanel()
    p.add_sequence(read_genbank(StringIO(SPLICED_CIRCULAR))[0])
    assert p._active is not None and len(p._active.annotations) == 3
    yield p
    p.deleteLater()


def _flip(panel) -> None:
    panel._op_combo.setCurrentIndex(panel._op_combo.findData("reverse_complement"))
    panel._apply_op()


def _saved(panel):
    """Whatever "Save all sequences as GenBank…" would write, re-read."""
    return SeqIO.read(StringIO(write_genbank([panel._active])), "genbank")


def _feature(record, ftype: str):
    return next(f for f in record.features if f.type == ftype)


def test_reverse_complement_through_the_panel_keeps_every_feature_s_bases(panel) -> None:
    before = SeqIO.read(StringIO(SPLICED_CIRCULAR), "genbank")
    want = [str(f.extract(before.seq)) for f in before.features]

    _flip(panel)
    saved = _saved(panel)
    wrong = [
        (i, f.type, str(f.location))
        for i, f in enumerate(saved.features)
        if str(f.extract(saved.seq)) != want[i]
    ]
    assert wrong == [], f"{len(wrong)} features read different bases after the flip: {wrong}"


def test_reverse_complement_through_the_panel_keeps_the_spliced_protein(panel) -> None:
    """One click, and the saved CDS still starts at its start codon."""
    _flip(panel)
    saved = _saved(panel)
    bases = str(_feature(saved, "CDS").extract(saved.seq))
    assert bases.startswith("ATG"), bases[:12]
    assert bases == SPLICED_CDS
    assert str(Seq(bases).translate()) == SPLICED_PROTEIN
    assert str(_feature(saved, "rep_origin").extract(saved.seq)) == ORIGIN_SPANNING


def test_the_feature_table_still_lists_every_feature_after_the_flip(panel) -> None:
    """A control: this passed with the defect present, so it is not what pins it."""
    _flip(panel)
    assert panel._feature_table.rowCount() == 3
    assert panel._active.length == 120
    assert panel._active.seq == str(
        Seq(str(read_genbank(StringIO(SPLICED_CIRCULAR))[0].seq)).reverse_complement()
    )
