"""CO-8, CO-9, CO-12 through the widget the user actually clicks.

The arithmetic tested in test_coords_edit_shift.py and
test_coords_reverse_complement.py is reached from exactly one place in the
shipped app: `Sequence` tab → Edit Operation → Apply
(``gui/panels/sequence_panel.py:_apply_op`` → ``_commit_edit``). This file
drives that widget headlessly — set the fields, press Apply — and then saves
the record the way "Save all sequences as GenBank…" does
(``gui/main_window.py:210`` calls ``write_genbank``), so the assertion is on
the file a user would end up with.

Oracle: the unedited record's extracted subsequences, from ``Bio.SeqIO`` on the
pinned bytes.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from conftest import bio_extract, requirement, spliced, touched  # noqa: E402

CO_8 = requirement("CO-8")
CO_12 = requirement("CO-12")
INSERT = CO_8["parameters"]["edits"][0]          # insert 30 bases at 3000
DELETE = CO_8["parameters"]["edits"][2]          # delete [2100, 2200)


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture
def panel(qapp, jcv):
    from gui.panels.sequence_panel import SequencePanel

    p = SequencePanel()
    p.add_sequence(jcv)
    assert p._active is not None and p._active.id == jcv.id
    assert len(p._active.annotations) == 28
    yield p
    p.deleteLater()


def _choose(panel, op: str) -> None:
    panel._op_combo.setCurrentIndex(panel._op_combo.findData(op))


def _screenshot(panel, name: str) -> None:
    """Grab the widget the way CLAUDE.md's GUI-verification recipe does."""
    out = os.environ.get("GX_SCREENSHOT_DIR")
    if out:
        panel.grab().save(os.path.join(out, f"{name}.png"))


def _saved_record(panel):
    """Whatever "Save all sequences as GenBank…" would write, re-read."""
    from backend.formats.genbank import write_genbank

    from conftest import reread

    return reread(write_genbank([panel._active]))


def test_co_8_insert_through_the_panel_leaves_untouched_features_alone(panel, jcv_bio):
    before = [bio_extract(f, jcv_bio) for f in jcv_bio.features]
    parts_before = [[tuple(p) for p in a.parts] for a in panel._active.annotations]

    _choose(panel, "insert")
    panel._pos_start.setValue(INSERT["start"])
    panel._insert_seq.setText(INSERT["inserted_bases"])
    panel._apply_op()
    _screenshot(panel, "co8-insert")

    assert len(panel._active.seq) == 5130 + INSERT["delta"]
    wrong = []
    for j, ann in enumerate(panel._active.annotations):
        if touched(parts_before[j], INSERT["start"], INSERT["end"]):
            continue
        got = spliced(panel._active.seq, [tuple(p) for p in ann.parts], ann.strand.value)
        if got != before[j]:
            wrong.append((j, ann.feature_type))
    assert wrong == [], wrong


def test_co_9_the_saved_file_after_an_insert_still_describes_the_same_features(panel, jcv_bio):
    before = [bio_extract(f, jcv_bio) for f in jcv_bio.features]
    parts_before = [[tuple(p) for p in a.parts] for a in panel._active.annotations]

    _choose(panel, "insert")
    panel._pos_start.setValue(INSERT["start"])
    panel._insert_seq.setText(INSERT["inserted_bases"])
    panel._apply_op()

    saved = _saved_record(panel)
    wrong = [
        (i, f.type, str(f.location))
        for i, f in enumerate(saved.features)
        if not touched(parts_before[i], INSERT["start"], INSERT["end"])
        and bio_extract(f, saved) != before[i]
    ]
    assert wrong == [], wrong


def test_co_9_delete_through_the_panel_shifts_the_downstream_features(panel, jcv_bio):
    before = [bio_extract(f, jcv_bio) for f in jcv_bio.features]
    parts_before = [[tuple(p) for p in a.parts] for a in panel._active.annotations]

    _choose(panel, "delete")
    panel._pos_start.setValue(DELETE["start"])
    panel._pos_end.setValue(DELETE["end"])
    panel._apply_op()
    _screenshot(panel, "co9-delete")

    assert len(panel._active.seq) == 5130 + DELETE["delta"]
    saved = _saved_record(panel)
    kept = [a for a in panel._active.annotations]
    assert len(saved.features) == len(kept)
    wrong = [
        (i, f.type)
        for i, f in enumerate(saved.features)
        if not touched(parts_before[i], DELETE["start"], DELETE["end"])
        and bio_extract(f, saved) != before[i]
    ]
    assert wrong == [], wrong


def test_co_12_reverse_complement_through_the_panel_keeps_every_feature_s_bases(panel, jcv_bio):
    before = [bio_extract(f, jcv_bio) for f in jcv_bio.features]

    _choose(panel, "reverse_complement")
    panel._apply_op()
    _screenshot(panel, "co12-reverse-complement")

    saved = _saved_record(panel)
    wrong = [
        (i, f.type, str(f.location))
        for i, f in enumerate(saved.features)
        if bio_extract(f, saved) != before[i]
    ]
    assert wrong == [], f"{len(wrong)} features read different bases after Reverse Complement: {wrong}"


def test_co_12_the_feature_table_still_lists_every_feature_after_the_flip(panel):
    """A user-visible check that the flip did not quietly lose annotations."""
    _choose(panel, "reverse_complement")
    panel._apply_op()
    assert panel._feature_table.rowCount() == 28
