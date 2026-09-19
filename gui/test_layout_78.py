"""#78 — the layout defects a design-review pass found, as geometry assertions.

Every panel in #78 passes its own tests; the defects are in where Qt put things,
not in what the code computed. So these tests build the *real* `MainWindow` at the
sizes the issue names, under the offscreen platform plugin, and assert on widget
geometry after the layout has actually run.

The acceptance list on the issue is the spec:

1. all ten tabs reachable without scrolling at 1500x950 and 1280x800
2. the icon reads as one unit with its label (a tab is a row, not a rotated column)
3. no CRISPR guide-table header is narrower than the text it renders
4. the TARGET SEQUENCE box does not hoard vertical space the results table needs
4b. the same for the Primers TEMPLATE SEQUENCE box at 1400x800
5. the ADD SEQUENCE MANUALLY title clears its own border, and the empty sequence
   list says it is empty

Items 2 and 5's border overlap are also judged by eye from the screenshots the
fixer saved; the assertions here are the part a machine can re-run.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import (  # noqa: E402
    QApplication, QStyle, QStyleOptionGroupBox, QToolButton,
)

from gui.main_window import MainWindow  # noqa: E402
from backend.models.sequence import MoleculeType, Sequence  # noqa: E402


#: The ten tabs, in the order the rail must keep showing them.
TAB_TEXT = [
    "🧬 Sequences", "🔍 DB Search", "✂️ CRISPR", "↔ Alignment", "🔩 Primers",
    "🔬 Variants", "🗺️ Pathways", "⛓ Workflow", "🧪 Structure", "⚙ Settings",
]

#: QHeaderView::section in gui/styles.py: padding 6px 8px, plus a 1px separator.
HEADER_PADDING = 16

DEMO_SEQ = ("ATGGCCTGTGGGCATCACGATGGCCTGTGGGAAACCTTTGGCAGATCCGTAGCTAGCTAGG" * 12)[:732]


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


def _window(app: QApplication, width: int, height: int, *, seeded: bool = True) -> MainWindow:
    """A shown MainWindow whose layout has run, with a sequence loaded."""
    # The startup update check opens a network connection from a QThread; no
    # layout test needs it, and it would outlive the window. Since #84 that is a
    # constructor argument, so this no longer has to monkeypatch the class.
    win = MainWindow(check_updates=False)
    win.resize(width, height)
    win.show()
    app.processEvents()
    if seeded:
        win._seq_panel.add_sequence(
            Sequence(id="TEST732", seq=DEMO_SEQ, molecule_type=MoleculeType.DNA,
                     description="732 nt demo")
        )
    app.processEvents()
    return win


def _show_tab(app: QApplication, win: MainWindow, panel) -> None:
    win._tabs.setCurrentWidget(panel)
    app.processEvents()


# ── 1. every tab reachable without scrolling ────────────────────────────────

@pytest.mark.parametrize("width,height", [(1500, 950), (1280, 800)])
def test_all_ten_tabs_fit_in_the_bar_without_scrolling(
    app: QApplication, width: int, height: int
) -> None:
    win = _window(app, width, height)
    bar = win._tabs.tabBar()

    assert bar.count() == 10, "all ten tabs must stay in the rail"

    outside = [
        (i, bar.tabText(i), bar.tabRect(i))
        for i in range(bar.count())
        if not bar.rect().contains(bar.tabRect(i))
    ]
    assert not outside, (
        f"at {width}x{height} these tabs fall outside the visible tab bar "
        f"{bar.rect()}: {outside} (#78.1)"
    )

    visible_scrollers = [
        b.objectName() for b in bar.findChildren(QToolButton) if b.isVisible()
    ]
    assert not visible_scrollers, (
        f"at {width}x{height} the rail shows scroll buttons {visible_scrollers}; "
        "a user must not have to find an arrow to reach Settings (#78.1)"
    )


# ── 2. icon and label read as one tab ───────────────────────────────────────

def test_a_tab_is_a_row_not_a_rotated_column(app: QApplication) -> None:
    """A tab wider than it is tall means upright text with its icon beside it.

    Qt's West-position default rotates the label 90 degrees, which stacks the icon
    and the text along the rail and makes ten tabs read as twenty items (#78.2).
    """
    win = _window(app, 1500, 950)
    bar = win._tabs.tabBar()
    metrics = bar.fontMetrics()

    for i in range(bar.count()):
        rect, text = bar.tabRect(i), bar.tabText(i)
        assert rect.width() >= metrics.horizontalAdvance(text), (
            f"tab {i} ({text!r}) is {rect.width()} px wide, too narrow for its "
            f"{metrics.horizontalAdvance(text)} px label (#78.2)"
        )
        # One line of text plus padding. A rotated label makes the tab as tall as
        # the label is long, which is where the ~100 px tabs came from.
        assert rect.height() <= metrics.height() + 20, (
            f"tab {i} ({text!r}) is {rect.height()} px tall for a "
            f"{metrics.height()} px line — its label runs down the rail instead of "
            "across it, which splits the icon off as a separate item (#78.2)"
        )


def test_tab_order_and_text_are_unchanged(app: QApplication) -> None:
    win = _window(app, 1500, 950)
    assert [win._tabs.tabText(i) for i in range(win._tabs.count())] == TAB_TEXT


# ── 3. the on-target score header is not truncated ──────────────────────────

def test_no_crispr_header_section_is_narrower_than_its_label(app: QApplication) -> None:
    win = _window(app, 1500, 950)
    _show_tab(app, win, win._crispr_panel)
    table = win._crispr_panel._table
    header = table.horizontalHeader()
    metrics = header.fontMetrics()

    too_narrow = []
    for col in range(table.columnCount()):
        # The stylesheet renders header sections upper-cased, so that is the
        # text whose width has to fit.
        label = table.horizontalHeaderItem(col).text().upper()
        needed = metrics.horizontalAdvance(label) + HEADER_PADDING
        if header.sectionSize(col) < needed:
            too_narrow.append((col, label, header.sectionSize(col), needed))

    assert not too_narrow, (
        "these guide-table headers are cut off (col, label, width, needed): "
        f"{too_narrow} (#78.3)"
    )


def _sections_narrower_than_labels(table) -> list[tuple[int, str, int, int]]:
    """The same criterion as above, for any table: (col, label, width, needed)."""
    header = table.horizontalHeader()
    metrics = header.fontMetrics()
    return [
        (col, table.horizontalHeaderItem(col).text().upper(), header.sectionSize(col),
         metrics.horizontalAdvance(table.horizontalHeaderItem(col).text().upper())
         + HEADER_PADDING)
        for col in range(table.columnCount())
        if header.sectionSize(col)
        < metrics.horizontalAdvance(table.horizontalHeaderItem(col).text().upper())
        + HEADER_PADDING
    ]


@pytest.mark.parametrize("ui_font_px", [15, 17, 19])
def test_result_headers_fit_their_labels_in_a_wider_ui_font(
    app: QApplication, ui_font_px: int
) -> None:
    """Item 3 on a platform whose UI font is not this container's.

    `ResizeToContents` alone sizes a section from the font Qt renders the
    *section* in, and from the label as the model stores it; the criterion above
    is the font `header.fontMetrics()` reports, upper-cased. Those two happen to
    agree within 3 px here and disagreed by 3 px the other way on the verifier's
    macOS run, which is how round 1 passed on Linux and still clipped ON-TARGET
    and OFF-TARGETS for a user. Widening the header's own font relative to the
    section's reproduces that platform difference deterministically, so the fix
    can be checked without a Mac.
    """
    win = _window(app, 1500, 950)
    for panel, table in ((win._crispr_panel, "_table"),
                         (win._primers_panel, "_pcr_table")):
        _show_tab(app, win, panel)
        widget = getattr(panel, table)
        widget.horizontalHeader().setStyleSheet(
            f"QHeaderView {{ font-size: {ui_font_px}px; }}"
        )
        app.processEvents()
        assert not _sections_narrower_than_labels(widget), (
            f"at a {ui_font_px}px UI font these {table} headers are cut off "
            f"(col, label, width, needed): "
            f"{_sections_narrower_than_labels(widget)} (#78.3)"
        )


# ── 4 / 4b. the input box must not hoard the results table's space ──────────

def test_crispr_target_box_keeps_to_its_size_hint(app: QApplication) -> None:
    win = _window(app, 1500, 950)
    _show_tab(app, win, win._crispr_panel)
    panel = win._crispr_panel
    box = panel._target_input.parentWidget()

    allowed = box.layout().sizeHint().height() + 40
    assert box.height() <= allowed, (
        f"TARGET SEQUENCE is {box.height()} px tall for a layout that asks for "
        f"{box.layout().sizeHint().height()} px (#78.4)"
    )
    assert panel._table.height() > 2 * box.height(), (
        f"the guide table ({panel._table.height()} px) should get the space the "
        f"target box ({box.height()} px) was holding (#78.4)"
    )


def test_primers_template_box_keeps_to_its_size_hint(app: QApplication) -> None:
    win = _window(app, 1400, 800)
    _show_tab(app, win, win._primers_panel)
    panel = win._primers_panel
    box = panel._template.parentWidget()

    allowed = box.layout().sizeHint().height() + 40
    assert box.height() <= allowed, (
        f"TEMPLATE SEQUENCE is {box.height()} px tall for a layout that asks for "
        f"{box.layout().sizeHint().height()} px (#78.4b)"
    )
    assert panel._pcr_table.height() > box.height(), (
        f"the primer table ({panel._pcr_table.height()} px) should get the space "
        f"the template box ({box.height()} px) was holding (#78.4b)"
    )


# ── 5. minor: group-box title, empty-list placeholder ───────────────────────

def test_group_box_title_does_not_overlap_its_border(app: QApplication) -> None:
    win = _window(app, 1500, 950)
    _show_tab(app, win, win._seq_panel)
    box = win._seq_panel._paste_id.parentWidget()
    assert box.title().lower() == "add sequence manually"

    option = QStyleOptionGroupBox()
    box.initStyleOption(option)
    style = box.style()
    label = style.subControlRect(
        QStyle.ComplexControl.CC_GroupBox, option, QStyle.SubControl.SC_GroupBoxLabel, box)
    frame = style.subControlRect(
        QStyle.ComplexControl.CC_GroupBox, option, QStyle.SubControl.SC_GroupBoxFrame, box)

    assert label.bottom() < frame.top(), (
        f"the {box.title()!r} title (bottom {label.bottom()}) runs into its own "
        f"border (top {frame.top()}) (#78.5)"
    )


def test_empty_sequence_list_says_it_is_empty(app: QApplication) -> None:
    win = _window(app, 1500, 950, seeded=False)
    _show_tab(app, win, win._seq_panel)
    panel = win._seq_panel

    assert panel._list.count() == 0
    placeholder = panel._list_placeholder
    assert placeholder.isVisible(), (
        "the empty sequence list shows no placeholder, so its 400 px of blank "
        "space reads as broken rather than empty (#78.5)"
    )
    assert placeholder.text().strip(), "the placeholder must actually say something"

    panel.add_sequence(
        Sequence(id="TEST732", seq=DEMO_SEQ, molecule_type=MoleculeType.DNA, description="")
    )
    app.processEvents()
    assert not placeholder.isVisible(), (
        "the placeholder must go away once a sequence is loaded"
    )
