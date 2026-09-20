"""A table header that is never narrower than the text it actually draws.

Round 1 of #78 fixed the truncated `ON-TARGET SCORE` header by switching the
guide table to `ResizeToContents`. That is still font-dependent, and it depended
on the wrong font: Qt sizes a section from the *section's* font and from the
string as the model stores it ("On-target"), while `gui/styles.py` draws header
sections upper-cased ("ON-TARGET") and the header widget's own font — the one
`header.fontMetrics()` reports — can be wider than the one Qt measured.

On this container that left 3 px of slack; on macOS it was 3 px short, which is
how round 2 came back. So do not trust either font: take whatever Qt computed and
raise it to the width the upper-cased label needs in the live `fontMetrics()`,
measured at the moment the section is sized. That holds on any platform, in any
font, at any DPI, without pinning a number.
"""

from __future__ import annotations

from PyQt6.QtCore import QEvent, QSize, Qt
from PyQt6.QtGui import QFontMetrics
from PyQt6.QtWidgets import QHeaderView, QStyle, QTableWidget


class LabelFittedHeaderView(QHeaderView):
    """Horizontal `QHeaderView` whose sections always fit their own labels."""

    #: `QHeaderView::section` in gui/styles.py is `padding: 6px 8px`, plus the
    #: 1 px separator between sections. The tests use the same number.
    SECTION_PADDING = 16

    #: The modes Qt sizes from a fixed number rather than from the content, so
    #: `sectionSizeFromContents` below never gets a say in how wide they are.
    _FIXED_WIDTH_MODES = (
        QHeaderView.ResizeMode.Interactive,
        QHeaderView.ResizeMode.Fixed,
    )

    def __init__(self, parent=None) -> None:
        super().__init__(Qt.Orientation.Horizontal, parent)
        self._enforcing = False
        #: Sections taken out of `Stretch` because the share it gave them was
        #: narrower than their label; put back as soon as there is room again.
        self._unstretched: set[int] = set()
        self.sectionCountChanged.connect(self._on_section_count_changed)

    # ── sizing ────────────────────────────────────────────────────────────
    def _label_width(self, logical_index: int) -> int:
        """Width the label needs in the fonts this header could render it in."""
        model = self.model()
        if model is None:
            return 0
        text = model.headerData(
            logical_index, self.orientation(), Qt.ItemDataRole.DisplayRole
        )
        if text is None:
            return 0
        text = str(text)

        # Both cases, because the stylesheet upper-cases what it draws; and both
        # fonts, because the section may be styled apart from the widget.
        candidates = [text, text.upper()]
        metrics = [self.fontMetrics()]
        font = model.headerData(
            logical_index, self.orientation(), Qt.ItemDataRole.FontRole
        )
        if font is not None:
            metrics.append(QFontMetrics(font))
        width = max(fm.horizontalAdvance(t) for fm in metrics for t in candidates)

        style = self.style()
        extra = self.SECTION_PADDING + 2 * style.pixelMetric(
            QStyle.PixelMetric.PM_HeaderMargin, None, self
        )
        if self.isSortIndicatorShown():
            extra += style.pixelMetric(QStyle.PixelMetric.PM_HeaderMarkSize, None, self)
        return width + extra

    def sectionSizeFromContents(self, logicalIndex: int) -> QSize:  # noqa: N803 (Qt)
        size = super().sectionSizeFromContents(logicalIndex)
        return QSize(max(size.width(), self._label_width(logicalIndex)), size.height())

    # ── sections Qt sizes from a number, not from their contents ──────────
    def _enforce_label_widths(self) -> None:
        """Widen — never narrow — the sections `sectionSizeFromContents` misses.

        `ResizeToContents` routes through the override above and `Stretch` takes
        whatever space is left, but an `Interactive` or `Fixed` section keeps
        `defaultSectionSize` (100 px) no matter how long its label is. That is
        how Structure's `RELATIVE SASA` — 116 px of text in a 100 px section —
        rendered as `ELATIVE SAS` (#98). Only ever grow, so no column ends up
        narrower than Qt would have made it.
        """
        if self._enforcing or self.model() is None:
            return
        self._enforcing = True
        try:
            for index in range(self.count()):
                if self.sectionResizeMode(index) not in self._FIXED_WIDTH_MODES:
                    continue
                needed = self._label_width(index)
                if self.sectionSize(index) < needed:
                    self.resizeSection(index, needed)
            self._rebalance_stretch()
        finally:
            self._enforcing = False

    # ── stretched sections, which get whatever space is left over ─────────
    def _stretch_share(self, sections: list[int]) -> int:
        """What each of `sections` would get if they were the stretched ones."""
        if not sections:
            return 0
        others = sum(
            self.sectionSize(i)
            for i in range(self.count())
            if i not in sections and not self.isSectionHidden(i)
        )
        return (self.width() - others) // len(sections)

    def _rebalance_stretch(self) -> None:
        """Take a starved column out of `Stretch`, and put it back when it fits.

        `Stretch` divides the space the other columns did not use, so it can hand
        a column less than its own label needs — and `resizeSection` cannot help,
        because the next layout recomputes the share. Qt offers no per-section
        floor for it (`minimumSectionSize` is header-wide and is not applied to
        stretched sections at all), so the only honest move is to stop stretching
        that column and give it the width its label needs, letting the table
        scroll. The two conditions below are exact inverses of each other, so a
        column cannot oscillate between them at a fixed width.
        """
        stretch = [
            i for i in range(self.count())
            if self.sectionResizeMode(i) == QHeaderView.ResizeMode.Stretch
        ]
        starved = [i for i in stretch if self.sectionSize(i) < self._label_width(i)]
        if starved:
            for index in stretch:
                self._unstretched.add(index)
                width = max(self.sectionSize(index), self._label_width(index),
                            self.defaultSectionSize())
                self.setSectionResizeMode(index, QHeaderView.ResizeMode.Interactive)
                self.resizeSection(index, width)
            return

        held = sorted(self._unstretched)
        if held and self._stretch_share(held) >= max(
            self._label_width(i) for i in held
        ):
            for index in held:
                self.setSectionResizeMode(index, QHeaderView.ResizeMode.Stretch)
            self._unstretched.clear()

    def resizeEvent(self, event) -> None:
        """A narrower table gives stretched sections a smaller share; re-check."""
        super().resizeEvent(event)
        self._enforce_label_widths()

    def _on_section_count_changed(self, _old: int, _new: int) -> None:
        self._enforce_label_widths()

    def _on_header_data_changed(self, orientation, _first: int, _last: int) -> None:
        if orientation == self.orientation():
            self._enforce_label_widths()

    def setModel(self, model) -> None:
        """Track header-label changes, so columns filled later are fitted too.

        Alignment builds its identity table's columns at results time, not at
        construction; without this it would be measured once, while empty (#98).
        """
        previous = self.model()
        if previous is not None:
            try:
                previous.headerDataChanged.disconnect(self._on_header_data_changed)
            except TypeError:  # never connected, e.g. a header swapped in twice
                pass
        super().setModel(model)
        if model is not None:
            model.headerDataChanged.connect(self._on_header_data_changed)
        self._enforce_label_widths()

    # ── keep the widths honest when the font underneath changes ───────────
    def changeEvent(self, event: QEvent) -> None:
        super().changeEvent(event)
        if event.type() in (
            QEvent.Type.FontChange,
            QEvent.Type.StyleChange,
            QEvent.Type.ApplicationFontChange,
        ):
            # Sections sized before a font/style swap were measured in the old
            # metrics; ask for them again rather than leave a stale width.
            self.resizeSections()
            self._enforce_label_widths()


def fit_header_to_labels(table, *, stretch: tuple[int, ...] = ()) -> LabelFittedHeaderView:
    """Give `table` a label-fitted header, sizing to contents but for `stretch`.

    Call it after the header labels are set, so the first size pass has text to
    measure.
    """
    header = table.horizontalHeader()
    if not isinstance(header, LabelFittedHeaderView):
        header = LabelFittedHeaderView(table)
        table.setHorizontalHeader(header)
    header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
    for col in stretch:
        header.setSectionResizeMode(col, QHeaderView.ResizeMode.Stretch)
    return header


def fitted_table(*args, **kwargs) -> QTableWidget:
    """A `QTableWidget` whose header fits its labels from construction onward.

    #78 installed `LabelFittedHeaderView` on two tables by hand and the other ten
    were missed, which is the whole of #98. Build results tables through here
    instead: the header is in place before any label or resize mode is set, so a
    panel gets the fitting without having to remember to ask for it.
    """
    table = QTableWidget(*args, **kwargs)
    table.setHorizontalHeader(LabelFittedHeaderView(table))
    return table
