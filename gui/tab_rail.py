"""The left tab rail: vertical stack, horizontal labels."""

from __future__ import annotations

from PyQt6.QtCore import QSize
from PyQt6.QtWidgets import QStyle, QStyleOptionTab, QStylePainter, QTabBar


class RailTabBar(QTabBar):
    """A West-position tab bar that draws each tab as one upright row.

    Qt's default for a West tab bar rotates the label 90 degrees, which makes a
    tab as *tall* as its label is long. Ten tabs then needed roughly 1100 px and
    the last two fell off a 950 px window behind a scroll arrow, so Settings was
    unreachable unless the user spotted it (#78.1). The rotation also laid the
    emoji out along the rail rather than beside the word, so the rail read as
    twenty small items instead of ten tabs (#78.2).

    Both come from the same place: the label's orientation. The shape is still
    drawn with the real (vertical) shape so the stylesheet's selected/hover rules
    apply, but the label is drawn as if the tab were at the top of the window,
    which is what upright means to the style. A tab is then one line tall, and
    ten of them fit any window this app supports.
    """

    #: Wide enough for the longest label, so the rail's width does not jump.
    MIN_WIDTH = 150
    _H_PADDING = 28
    _V_PADDING = 16

    def tabSizeHint(self, index: int) -> QSize:
        metrics = self.fontMetrics()
        width = metrics.horizontalAdvance(self.tabText(index)) + self._H_PADDING
        return QSize(max(width, self.MIN_WIDTH), metrics.height() + self._V_PADDING)

    def paintEvent(self, event) -> None:  # noqa: ARG002 — Qt signature
        painter = QStylePainter(self)
        option = QStyleOptionTab()
        for index in range(self.count()):
            self.initStyleOption(option, index)
            painter.drawControl(QStyle.ControlElement.CE_TabBarTabShape, option)
            # RoundedNorth is the "label runs left to right" case; the rect stays
            # the real one, so only the text and icon orientation changes.
            option.shape = QTabBar.Shape.RoundedNorth
            painter.drawControl(QStyle.ControlElement.CE_TabBarTabLabel, option)
