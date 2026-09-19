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
from PyQt6.QtWidgets import QHeaderView, QStyle


class LabelFittedHeaderView(QHeaderView):
    """Horizontal `QHeaderView` whose sections always fit their own labels."""

    #: `QHeaderView::section` in gui/styles.py is `padding: 6px 8px`, plus the
    #: 1 px separator between sections. The tests use the same number.
    SECTION_PADDING = 16

    def __init__(self, parent=None) -> None:
        super().__init__(Qt.Orientation.Horizontal, parent)

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


def fit_header_to_labels(table, *, stretch: tuple[int, ...] = ()) -> LabelFittedHeaderView:
    """Give `table` a label-fitted header, sizing to contents but for `stretch`.

    Call it after the header labels are set, so the first size pass has text to
    measure.
    """
    header = LabelFittedHeaderView(table)
    table.setHorizontalHeader(header)
    header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
    for col in stretch:
        header.setSectionResizeMode(col, QHeaderView.ResizeMode.Stretch)
    return header
