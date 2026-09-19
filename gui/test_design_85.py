"""#85 — three design defects found by eye at 1280x800 and 1500x950, as assertions.

Same method as `test_layout_78.py`: build the *real* `MainWindow` under the
offscreen platform plugin at the sizes the issue names, let the layout actually
run, and assert on what Qt did. The issue's acceptance list is the spec:

1. at 1280x800 every label in the PCR **Parameters** row renders its full text —
   its width is at least `fontMetrics().horizontalAdvance(text)`, measured from
   the font that actually draws *that* label (#78 round 2's lesson: a width that
   fits this container's font clipped on the verifier's macOS font)
2. `Remove` (Sequences) is disabled with nothing selected and enabled with a
   selection; `Delete` (Primers preset) is disabled while no deletable saved
   preset is selected
3. the ADD SEQUENCE MANUALLY rows sit at their size-hint spacing, with no
   stretch opened up between them
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import json  # noqa: E402

import pytest  # noqa: E402
from PyQt6.QtGui import QFont  # noqa: E402
from PyQt6.QtWidgets import (  # noqa: E402
    QApplication, QFormLayout, QLabel, QPushButton, QWidget,
)

from gui.main_window import MainWindow  # noqa: E402
from gui.panels import primers_panel  # noqa: E402
from backend.models.sequence import MoleculeType, Sequence  # noqa: E402

DEMO_SEQ = ("ATGGCCTGTGGGCATCACGATGGCCTGTGGGAAACCTTTGGCAGATCCGTAGCTAGCTAGG" * 12)[:732]

#: What the placeholder row of the preset combo says when nothing is chosen.
NO_PRESET = "— select preset —"


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


def _window(app: QApplication, width: int, height: int) -> MainWindow:
    """A shown MainWindow whose layout has run, with an empty sequence list.

    `check_updates=False` since #84 — no layout test wants the startup update
    check's QThread outliving the window.
    """
    win = MainWindow(check_updates=False)
    win.resize(width, height)
    win.show()
    app.processEvents()
    return win


def _show_tab(app: QApplication, win: MainWindow, panel: QWidget) -> None:
    win._tabs.setCurrentWidget(panel)
    app.processEvents()


def _button(panel: QWidget, text: str) -> QPushButton:
    """The panel's button with this exact label, found the way a user finds it."""
    matches = [b for b in panel.findChildren(QPushButton) if b.text() == text]
    assert len(matches) == 1, f"expected exactly one {text!r} button, found {len(matches)}"
    return matches[0]


def _demo_sequence(seq_id: str = "TEST732") -> Sequence:
    return Sequence(id=seq_id, seq=DEMO_SEQ, molecule_type=MoleculeType.DNA,
                    description="732 nt demo")


# ── 1. the PCR Parameters labels render their full text ─────────────────────

def _clipped_labels(group: QWidget) -> list[tuple[str, int, int]]:
    """(text, width, needed) for every label in `group` narrower than its text.

    `label.fontMetrics()` is the font Qt will actually draw this label with,
    stylesheet and inherited font included — not the panel's or the app's.
    """
    return [
        (lb.text(), lb.width(), lb.fontMetrics().horizontalAdvance(lb.text()))
        for lb in group.findChildren(QLabel)
        if lb.text() and lb.width() < lb.fontMetrics().horizontalAdvance(lb.text())
    ]


def _params_group(win: MainWindow) -> QWidget:
    group = win._primers_panel._prod_min.parentWidget()
    assert group.title().lower() == "parameters"
    return group


@pytest.mark.parametrize("width,height", [(1280, 800), (1500, 950)])
def test_every_pcr_parameter_label_renders_its_full_text(
    app: QApplication, width: int, height: int
) -> None:
    """#85.1 — 'Product min (bp):' lost its '):' to five other labels' widths."""
    win = _window(app, width, height)
    _show_tab(app, win, win._primers_panel)
    group = _params_group(win)

    assert not _clipped_labels(group), (
        f"at {width}x{height} these PCR Parameters labels are cut off "
        f"(text, width, needed): {_clipped_labels(group)} (#85.1)"
    )


@pytest.mark.parametrize("ui_font_pt", [11, 13, 15])
def test_pcr_parameter_labels_fit_in_a_wider_ui_font(
    app: QApplication, ui_font_pt: int
) -> None:
    """#85.1 on a platform whose UI font is not this container's.

    Round 2 of #78 passed on Linux and still clipped on macOS, because the fix
    was sized against one font. Widening the font the labels are *drawn* with
    reproduces that difference deterministically: whatever the layout does, it
    has to do it for the metrics of the font in front of it.
    """
    win = _window(app, 1280, 800)
    _show_tab(app, win, win._primers_panel)
    group = _params_group(win)

    font = QFont(group.font())
    font.setPointSize(ui_font_pt)
    group.setFont(font)
    group.layout().invalidate()
    group.layout().activate()
    app.processEvents()

    assert not _clipped_labels(group), (
        f"at a {ui_font_pt}pt UI font these PCR Parameters labels are cut off "
        f"(text, width, needed): {_clipped_labels(group)} (#85.1)"
    )


def test_the_six_pcr_parameter_spin_boxes_are_all_still_there(app: QApplication) -> None:
    """Nothing gets dropped to make the row fit."""
    win = _window(app, 1280, 800)
    _show_tab(app, win, win._primers_panel)
    panel = win._primers_panel

    for name in ("_prod_min", "_prod_max", "_len_min", "_len_max", "_tm_min", "_tm_max"):
        spin = getattr(panel, name)
        assert spin.isVisible(), f"{name} vanished from the Parameters row (#85.1)"


# ── 2. destructive actions are disabled until they have a target ────────────

def test_remove_is_disabled_until_a_sequence_is_selected(app: QApplication) -> None:
    """#85.2 — full-width bright red 'Remove' on an empty list can do nothing."""
    win = _window(app, 1500, 950)
    _show_tab(app, win, win._seq_panel)
    panel = win._seq_panel
    remove = _button(panel, "Remove")

    assert panel._list.count() == 0
    assert not remove.isEnabled(), (
        "'Remove' is enabled while the sequence list is empty — the most "
        "prominent control on the panel, with nothing to act on (#85.2)"
    )

    panel.add_sequence(_demo_sequence())
    app.processEvents()
    assert panel._list.currentItem() is not None
    assert remove.isEnabled(), (
        "'Remove' must come back once a sequence is selected (#85.2)"
    )

    remove.click()
    app.processEvents()
    assert panel._list.count() == 0
    assert not remove.isEnabled(), (
        "'Remove' stayed enabled after removing the last sequence (#85.2)"
    )


def test_remove_follows_the_selection_across_several_sequences(app: QApplication) -> None:
    """Enabled state tracks the selection, not just the list being non-empty."""
    win = _window(app, 1500, 950)
    _show_tab(app, win, win._seq_panel)
    panel = win._seq_panel
    remove = _button(panel, "Remove")

    panel.add_sequence(_demo_sequence("A"))
    panel.add_sequence(_demo_sequence("B"))
    app.processEvents()
    assert remove.isEnabled()

    panel._list.setCurrentItem(None)
    app.processEvents()
    assert not remove.isEnabled(), (
        "'Remove' is enabled with two sequences loaded and neither selected — "
        "there is still no target (#85.2)"
    )

    panel._list.setCurrentRow(1)
    app.processEvents()
    assert remove.isEnabled()


def test_remove_keeps_its_danger_styling(app: QApplication) -> None:
    """Disabling it is the fix; un-styling it is not (the issue says so)."""
    win = _window(app, 1500, 950)
    _show_tab(app, win, win._seq_panel)
    assert _button(win._seq_panel, "Remove").objectName() == "danger"


def test_delete_preset_is_disabled_while_no_saved_preset_is_selected(
    app: QApplication, monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    """#85.2 — 'Delete' offered against '— select preset —' and against built-ins.

    Neither is a target: the placeholder is not a preset, and deleting a
    built-in has always been refused with a warning box.
    """
    monkeypatch.setattr(primers_panel, "_PRESETS_FILE", tmp_path / "primer_presets.json")
    win = _window(app, 1500, 950)
    _show_tab(app, win, win._primers_panel)
    panel = win._primers_panel
    delete = _button(panel, "Delete")
    combo = panel._preset_combo

    panel._refresh_presets()
    app.processEvents()
    assert combo.currentText() == NO_PRESET
    assert not delete.isEnabled(), (
        f"'Delete' is enabled while the combo shows {NO_PRESET!r} (#85.2)"
    )

    combo.setCurrentText("RT-qPCR")
    app.processEvents()
    assert not delete.isEnabled(), (
        "'Delete' is enabled for the built-in preset 'RT-qPCR', which the panel "
        "refuses to delete anyway (#85.2)"
    )

    (tmp_path / "primer_presets.json").write_text(json.dumps(
        {"My bench mix": {"prod_min": 120, "prod_max": 900, "len_min": 19,
                          "len_max": 23, "tm_min": 57.0, "tm_max": 63.0}}))
    panel._refresh_presets()
    combo.setCurrentText("My bench mix")
    app.processEvents()
    assert delete.isEnabled(), (
        "'Delete' must be available for a preset the user saved themselves (#85.2)"
    )


# ── 3. the Add-sequence form sits at its natural height ─────────────────────

def _form(win: MainWindow) -> tuple[QWidget, QFormLayout]:
    box = win._seq_panel._paste_id.parentWidget()
    assert box.title().lower() == "add sequence manually"
    return box, box.layout()


@pytest.mark.parametrize("width,height", [(1500, 950), (1280, 800)])
def test_add_sequence_form_rows_sit_at_size_hint_spacing(
    app: QApplication, width: int, height: int
) -> None:
    """#85.3 — ID / Type / Seq spread over ~300 px of box with gaps between them."""
    win = _window(app, width, height)
    _show_tab(app, win, win._seq_panel)
    box, form = _form(win)

    fields = [
        form.itemAt(row, QFormLayout.ItemRole.FieldRole).widget()
        for row in range(form.rowCount())
        if form.itemAt(row, QFormLayout.ItemRole.FieldRole) is not None
        and form.itemAt(row, QFormLayout.ItemRole.FieldRole).widget() is not None
    ]
    assert len(fields) >= 4, "ID / Type / Seq / Add must all still be in the form"

    spacing = form.verticalSpacing()
    stretched = [
        (a.objectName() or type(a).__name__, b.objectName() or type(b).__name__,
         b.y() - (a.y() + a.height()))
        for a, b in zip(fields, fields[1:])
        if b.y() - (a.y() + a.height()) > spacing
    ]
    assert not stretched, (
        f"at {width}x{height} the form has stretch between rows, not its "
        f"{spacing} px spacing (previous, next, gap): {stretched} (#85.3)"
    )


def test_add_sequence_box_keeps_to_its_size_hint(app: QApplication) -> None:
    """The box should sit at its natural height, leaving the rest to the list."""
    win = _window(app, 1500, 950)
    _show_tab(app, win, win._seq_panel)
    box, _ = _form(win)

    assert box.height() <= box.sizeHint().height() + 2, (
        f"ADD SEQUENCE MANUALLY is {box.height()} px tall for a box that asks "
        f"for {box.sizeHint().height()} px, so the form stretches to fill it "
        "instead of sitting at its natural height (#85.3)"
    )


def test_the_sequence_list_gets_the_space_the_form_was_holding(app: QApplication) -> None:
    win = _window(app, 1500, 950)
    _show_tab(app, win, win._seq_panel)
    panel = win._seq_panel
    box, _ = _form(win)

    assert panel._list.height() > box.height(), (
        f"the sequence list ({panel._list.height()} px) should get the space the "
        f"add-sequence box ({box.height()} px) was holding (#85.3)"
    )
