"""#98 — every result table's header must fit the text it draws, not just two.

#78 round 2 added `LabelFittedHeaderView` and installed it by hand on the CRISPR
guide table and the Primers PCR table. Every other `QTableWidget` in `gui/` kept
the default header, so none of them measures its own label: Structure's
`RELATIVE SASA` renders as `ELATIVE SAS` at 1280x800 today, and the rest are
simply unmeasured — which of them clip depends on the platform's fonts, exactly
the way round 1 passed on Linux and failed on macOS.

So this walks *every* table in a built `MainWindow` rather than a hand-kept list,
at both sizes the issue names and at several UI font sizes, and applies the same
inequality #78 uses: a section is at least as wide as its upper-cased label in
`header.fontMetrics()` plus the stylesheet's padding.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication, QHeaderView, QTableWidget  # noqa: E402

from backend.models.sequence import MoleculeType, Sequence  # noqa: E402
from gui.main_window import MainWindow  # noqa: E402
from gui.table_header import LabelFittedHeaderView  # noqa: E402

#: `QHeaderView::section` in gui/styles.py is `padding: 6px 8px`, plus the 1 px
#: separator between sections — the same constant gui/test_layout_78.py uses.
HEADER_PADDING = 16

#: The sizes the issue names. A header that fits at one and clips at the other is
#: still a bug, so both are checked for every table.
SIZES = [(1280, 800), (1500, 950)]

#: More than one, so a width that happens to fit this container's metrics cannot
#: carry the test (#78 round 2's lesson).
UI_FONT_PX = [15, 17, 19]

DEMO_SEQ = ("ATGGCCTGTGGGCATCACGATGGCCTGTGGGAAACCTTTGGCAGATCCGTAGCTAGCTAGG" * 12)[:732]


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


def _window(app: QApplication, width: int, height: int) -> MainWindow:
    """A shown MainWindow whose layout has run, with a sequence loaded."""
    win = MainWindow(check_updates=False)
    win.resize(width, height)
    win.show()
    app.processEvents()
    win._seq_panel.add_sequence(
        Sequence(id="TEST732", seq=DEMO_SEQ, molecule_type=MoleculeType.DNA,
                 description="732 nt demo")
    )
    app.processEvents()
    return win


def _tables(win: MainWindow) -> list[QTableWidget]:
    """Every `QTableWidget` the window builds — discovered, never hand-listed.

    A hand-kept list is how #78 shipped the fitter on two tables and missed the
    other ten, so the walk is what keeps this honest as panels are added.
    """
    tables = win.findChildren(QTableWidget)
    assert tables, "no tables found in the window — the walk is measuring nothing"
    return tables


def _name(win: MainWindow, table: QTableWidget) -> str:
    """`PanelClass._attr` for a table, so a failure says which one broke."""
    parent = table.parentWidget()
    while parent is not None:
        for attr, value in vars(parent).items():
            if value is table:
                return f"{type(parent).__name__}{attr}"
        parent = parent.parentWidget()
    return repr(table)


def _show_every_tab(app: QApplication, win: MainWindow) -> None:
    """Lay out every tab, since a tab that was never current never laid out."""
    for index in range(win._tabs.count()):
        win._tabs.setCurrentIndex(index)
        app.processEvents()


def _clipped(table: QTableWidget) -> list[tuple[int, str, int, int]]:
    """Sections narrower than their own label: (col, label, width, needed)."""
    header = table.horizontalHeader()
    metrics = header.fontMetrics()
    out = []
    for col in range(table.columnCount()):
        item = table.horizontalHeaderItem(col)
        if item is None:
            continue
        # The stylesheet draws header sections upper-cased, so that is the text
        # whose width has to fit.
        label = item.text().upper()
        needed = metrics.horizontalAdvance(label) + HEADER_PADDING
        if header.sectionSize(col) < needed:
            out.append((col, label, header.sectionSize(col), needed))
    return out


# ── the acceptance walk ─────────────────────────────────────────────────────

@pytest.mark.parametrize("width,height", SIZES)
def test_every_table_header_fits_its_labels(
    app: QApplication, width: int, height: int
) -> None:
    """Structure's `RELATIVE SASA` is the one that clips today; none may."""
    win = _window(app, width, height)
    _show_every_tab(app, win)

    broken = {
        _name(win, table): _clipped(table)
        for table in _tables(win)
        if _clipped(table)
    }
    assert not broken, (
        f"at {width}x{height} these header sections are cut off "
        f"(col, label, width, needed): {broken} (#98)"
    )


@pytest.mark.parametrize("ui_font_px", UI_FONT_PX)
def test_every_table_header_fits_its_labels_in_a_wider_ui_font(
    app: QApplication, ui_font_px: int
) -> None:
    """The same walk on a platform whose UI font is not this container's.

    Widening the header widget's own font relative to the font Qt sized the
    section in reproduces the Linux/macOS metric gap deterministically, so a
    width that only fits one platform cannot pass (#78 round 2).
    """
    win = _window(app, 1280, 800)
    _show_every_tab(app, win)

    broken = {}
    for table in _tables(win):
        table.horizontalHeader().setStyleSheet(
            f"QHeaderView {{ font-size: {ui_font_px}px; }}"
        )
        app.processEvents()
        clipped = _clipped(table)
        if clipped:
            broken[_name(win, table)] = clipped

    assert not broken, (
        f"at a {ui_font_px}px UI font these header sections are cut off "
        f"(col, label, width, needed): {broken} (#98)"
    )


# ── the fitter is installed everywhere, including on tables filled later ────

def test_every_table_uses_the_label_fitted_header(app: QApplication) -> None:
    """Including Alignment's identity table, which has no columns until results.

    A table whose columns arrive later cannot be measured by the walk above, so
    the structural check is what keeps it from being silently unmeasured.
    """
    win = _window(app, 1280, 800)
    _show_every_tab(app, win)

    unfitted = [
        _name(win, table)
        for table in _tables(win)
        if not isinstance(table.horizontalHeader(), LabelFittedHeaderView)
    ]
    assert not unfitted, (
        f"these tables still use the default header, so nothing measures their "
        f"labels: {unfitted} (#98)"
    )


def test_headers_stay_fitted_when_columns_arrive_after_construction(
    app: QApplication,
) -> None:
    """Alignment fills its identity table at results time; it must fit then too."""
    win = _window(app, 1280, 800)
    _show_every_tab(app, win)

    table = win._align_panel._identity_table
    ids = ["NM_007294.4", "ENST00000357654", "XM_006712345.3"]
    table.setColumnCount(len(ids))
    table.setHorizontalHeaderLabels(ids)
    app.processEvents()

    assert not _clipped(table), (
        f"the identity table clips the accessions it is given "
        f"(col, label, width, needed): {_clipped(table)} (#98)"
    )


# ── nothing got narrower, lost, or reordered ────────────────────────────────

#: (attribute path on MainWindow, header labels) as the panels build them today.
#: Fitting a header must not drop, rename or reorder a column.
EXPECTED_COLUMNS = {
    "_seq_panel._feature_table":
        ["Feature Type", "Strand", "Start", "End", "Qualifiers"],
    "_pathways_panel._reactome_table":
        ["Pathway", "Entities Found", "Entities Total", "p-value", "FDR"],
    "_pathways_panel._string_table":
        ["Protein A", "Protein B", "Score", "Interaction Type"],
    "_variants_panel._table":
        ["CHROM", "POS", "REF", "ALT", "Type", "Gene", "Clinical Sig.",
         "gnomAD AF"],
    "_search_panel._table":
        ["Name", "Database", "Description", "Organism", "Accession"],
    "_crispr_panel._table": None,
    "_primers_panel._pcr_table": None,
    "_primers_panel._re_table": ["Enzyme", "Recognition", "Count", "Positions"],
    "_primers_panel._dig_table": ["Fragment", "Size (bp)", "Start", "End"],
    "_workflow_panel._steps_table": ["Step", "Parameters (JSON)"],
    "_workflow_panel._results_table": ["Step", "Status", "Detail"],
    "_structure_panel._table":
        ["Residue #", "Residue", "Chain", "Classification", "Relative SASA"],
}


def _resolve(win: MainWindow, path: str):
    obj = win
    for part in path.split("."):
        obj = getattr(obj, part)
    return obj


@pytest.mark.parametrize("path,labels", sorted(EXPECTED_COLUMNS.items()))
def test_no_table_loses_a_column_or_changes_order(
    app: QApplication, path: str, labels: list[str] | None
) -> None:
    win = _window(app, 1280, 800)
    _show_every_tab(app, win)
    table = _resolve(win, path)

    if labels is None:  # column set not pinned here; only that it is non-empty
        assert table.columnCount() > 0, f"{path} lost its columns (#98)"
        return

    actual = [
        table.horizontalHeaderItem(col).text() for col in range(table.columnCount())
    ]
    assert actual == labels, f"{path} changed its columns (#98)"


def test_no_section_is_narrower_than_the_default_header_would_make_it(
    app: QApplication,
) -> None:
    """The fitter only ever widens: no column ends up below Qt's own sizing.

    Qt's default is `defaultSectionSize` for an Interactive section and
    `sectionSizeHint` for a dynamic one; fitting raises those to the label width
    and must never lower them.
    """
    win = _window(app, 1280, 800)
    _show_every_tab(app, win)

    shrunk = []
    for table in _tables(win):
        header = table.horizontalHeader()
        for col in range(table.columnCount()):
            mode = header.sectionResizeMode(col)
            if mode == QHeaderView.ResizeMode.Stretch:
                continue  # stretch takes whatever is left; nothing to compare to
            floor = (
                header.defaultSectionSize()
                if mode in (QHeaderView.ResizeMode.Interactive,
                            QHeaderView.ResizeMode.Fixed)
                else header.sectionSizeHint(col)
            )
            if header.sectionSize(col) < floor:
                shrunk.append((_name(win, table), col, header.sectionSize(col), floor))

    assert not shrunk, (
        f"these sections came out narrower than Qt's own sizing "
        f"(table, col, width, floor): {shrunk} (#98)"
    )
