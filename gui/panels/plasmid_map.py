"""Circular plasmid map widget (epic #15).

This is the project's first QPainter-based custom widget — kept deliberately
simple and well-commented as the pattern later sub-issues (feature arcs #24,
cut-site ticks #25, click sync #26) and epic #17's assembly view will follow.

This skeleton (#23) draws only the backbone: an outer ring, bp-position tick
marks at round intervals, and the sequence name + size in the centre. Feature
and enzyme rendering are layered on in the follow-up sub-issues.

Coordinate convention: bp 0 sits at the top (12 o'clock) and position
increases clockwise, matching SnapGene/Benchling. For a bp position ``p`` the
angle is ``theta = 2*pi * p / length`` and the point on a circle of radius
``r`` is ``(cx + r*sin(theta), cy - r*cos(theta))``.
"""

from __future__ import annotations

import math

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QColor, QFont, QPainter, QPen
from PyQt6.QtWidgets import QWidget

from backend.models.sequence import Sequence

# Candidate tick spacings (bp); the smallest one giving <= ~12 ticks is used.
_TICK_STEPS = [50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000]


def _tick_interval(length: int) -> int:
    """Pick a round bp interval that yields a readable number of ticks."""
    for step in _TICK_STEPS:
        if length / step <= 12:
            return step
    return _TICK_STEPS[-1]


class PlasmidMapWidget(QWidget):
    """Custom-painted circular map of a plasmid backbone."""

    def __init__(self) -> None:
        super().__init__()
        self._name = ""
        self._length = 0
        self._is_circular = False
        self.setMinimumSize(240, 240)

    def set_sequence(self, seq: Sequence) -> None:
        """Point the map at a sequence (intended for circular ones)."""
        self._name = seq.name or seq.id
        self._length = seq.length
        self._is_circular = seq.is_circular
        self.update()

    def _point(self, cx: float, cy: float, radius: float, pos: int) -> QPointF:
        theta = 2 * math.pi * (pos / self._length) if self._length else 0.0
        return QPointF(cx + radius * math.sin(theta), cy - radius * math.cos(theta))

    def paintEvent(self, event) -> None:  # noqa: N802 (Qt override)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        w, h = self.width(), self.height()
        cx, cy = w / 2, h / 2
        radius = min(w, h) / 2 - 34  # leave a margin for tick labels

        if self._length <= 0 or radius <= 0:
            painter.end()
            return

        # ── Backbone ring ─────────────────────────────────────────────────
        painter.setPen(QPen(QColor("#94a3b8"), 2))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawEllipse(QPointF(cx, cy), radius, radius)

        # ── bp-position ticks + labels ────────────────────────────────────
        interval = _tick_interval(self._length)
        tick_pen = QPen(QColor("#64748b"), 1)
        label_font = QFont("JetBrains Mono,Fira Code,Courier New", 8)
        painter.setFont(label_font)
        for pos in range(0, self._length, interval):
            inner = self._point(cx, cy, radius - 6, pos)
            outer = self._point(cx, cy, radius + 6, pos)
            painter.setPen(tick_pen)
            painter.drawLine(inner, outer)
            # Label just outside the tick.
            label_pt = self._point(cx, cy, radius + 20, pos)
            painter.setPen(QPen(QColor("#cbd5e1"), 1))
            rect = QRectF(label_pt.x() - 24, label_pt.y() - 8, 48, 16)
            painter.drawText(rect, Qt.AlignmentFlag.AlignCenter, self._fmt_bp(pos))

        # ── Centre label: name + size + topology ──────────────────────────
        painter.setPen(QPen(QColor("#e2e8f0"), 1))
        painter.setFont(QFont("Inter,Helvetica,Arial", 11, QFont.Weight.Bold))
        name_rect = QRectF(cx - radius, cy - 16, 2 * radius, 18)
        painter.drawText(name_rect, Qt.AlignmentFlag.AlignCenter, self._name or "(unnamed)")
        painter.setPen(QPen(QColor("#94a3b8"), 1))
        painter.setFont(QFont("Inter,Helvetica,Arial", 9))
        topo = "circular" if self._is_circular else "linear"
        size_rect = QRectF(cx - radius, cy + 4, 2 * radius, 16)
        painter.drawText(size_rect, Qt.AlignmentFlag.AlignCenter, f"{self._length:,} bp · {topo}")
        painter.end()

    @staticmethod
    def _fmt_bp(pos: int) -> str:
        if pos == 0:
            return "1"
        if pos >= 1000 and pos % 1000 == 0:
            return f"{pos // 1000}k"
        return str(pos)
