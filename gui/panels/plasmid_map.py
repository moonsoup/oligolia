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
from PyQt6.QtGui import QBrush, QColor, QFont, QPainter, QPen, QPolygonF
from PyQt6.QtWidgets import QWidget

from backend.models.sequence import Annotation, Sequence, Strand
from gui.panels.feature_colors import feature_color_map

# Candidate tick spacings (bp); the smallest one giving <= ~12 ticks is used.
_TICK_STEPS = [50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000]


def _tick_interval(length: int) -> int:
    """Pick a round bp interval that yields a readable number of ticks."""
    for step in _TICK_STEPS:
        if length / step <= 12:
            return step
    return _TICK_STEPS[-1]


def _assign_lanes(annotations: list[Annotation]) -> list[tuple[Annotation, int]]:
    """Greedily pack features into concentric lanes so none overlap in a lane.

    Returns (annotation, lane_index) pairs; lane 0 is the outermost.
    """
    lane_ends: list[int] = []  # last end position occupied in each lane
    placed: list[tuple[Annotation, int]] = []
    for ann in sorted(annotations, key=lambda a: (a.start, a.end)):
        lane = next((i for i, end in enumerate(lane_ends) if ann.start >= end), None)
        if lane is None:
            lane = len(lane_ends)
            lane_ends.append(ann.end)
        else:
            lane_ends[lane] = ann.end
        placed.append((ann, lane))
    return placed


class PlasmidMapWidget(QWidget):
    """Custom-painted circular map of a plasmid backbone."""

    _LANE_WIDTH = 9  # arc thickness / lane spacing (px)

    def __init__(self) -> None:
        super().__init__()
        self._name = ""
        self._length = 0
        self._is_circular = False
        self._lanes: list[tuple[Annotation, int]] = []
        self._color_map: dict[str, QColor] = {}
        self.setMinimumSize(240, 240)

    def set_sequence(self, seq: Sequence) -> None:
        """Point the map at a sequence (intended for circular ones)."""
        self._name = seq.name or seq.id
        self._length = seq.length
        self._is_circular = seq.is_circular
        self._lanes = _assign_lanes(seq.annotations)
        self._color_map = feature_color_map(seq.annotations)
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

        # ── Feature arcs (issue #24) ──────────────────────────────────────
        self._draw_features(painter, cx, cy, radius)

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

    def _draw_features(self, painter: QPainter, cx: float, cy: float, radius: float) -> None:
        """Draw each annotation as a colored arc + strand arrowhead + label."""
        label_font = QFont("Inter,Helvetica,Arial", 8)
        for ann, lane in self._lanes:
            r = radius - 12 - lane * (self._LANE_WIDTH + 2)
            if r < radius * 0.35:  # too many overlapping lanes to fit; skip inner ones
                continue
            color = self._color_map.get(ann.feature_type, QColor("#94a3b8"))

            # Arc: Qt angles are 1/16 deg, 0 at 3 o'clock, CCW positive. Our bp 0
            # is at 12 o'clock increasing clockwise, so map bp p -> 90 - 360*p/L.
            span_bp = max(ann.end - ann.start, 1)
            start_angle = (90.0 - 360.0 * ann.start / self._length) * 16
            span_angle = -(360.0 * span_bp / self._length) * 16
            rect = QRectF(cx - r, cy - r, 2 * r, 2 * r)
            painter.setPen(QPen(color, self._LANE_WIDTH, Qt.PenStyle.SolidLine,
                                Qt.PenCapStyle.FlatCap))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawArc(rect, round(start_angle), round(span_angle))

            # Strand arrowhead: PLUS at the end (clockwise), MINUS at the start
            # (counter-clockwise); none for BOTH/unstranded.
            if ann.strand == Strand.PLUS:
                self._draw_arrowhead(painter, cx, cy, r, ann.end, +1, color)
            elif ann.strand == Strand.MINUS:
                self._draw_arrowhead(painter, cx, cy, r, ann.start, -1, color)

            # Label with a short leader line out to the arc's mid-angle.
            mid = (ann.start + ann.end) / 2
            label = (ann.qualifiers.get("label") or ann.qualifiers.get("gene")
                     or ann.qualifiers.get("product") or ann.feature_type)
            painter.setPen(QPen(QColor("#475569"), 1))
            painter.drawLine(self._point(cx, cy, r, mid),
                             self._point(cx, cy, radius + 12, mid))
            anchor = self._point(cx, cy, radius + 16, mid)
            painter.setFont(label_font)
            painter.setPen(QPen(color.lighter(130), 1))
            on_left = anchor.x() < cx
            box = QRectF(anchor.x() - (120 if on_left else 0), anchor.y() - 8, 120, 16)
            painter.drawText(
                box,
                (Qt.AlignmentFlag.AlignRight if on_left else Qt.AlignmentFlag.AlignLeft)
                | Qt.AlignmentFlag.AlignVCenter,
                str(label)[:24],
            )

    def _draw_arrowhead(self, painter: QPainter, cx: float, cy: float, r: float,
                        pos: int, direction: int, color: QColor) -> None:
        theta = 2 * math.pi * (pos / self._length) if self._length else 0.0
        # Unit tangent for increasing bp (clockwise), then oriented by strand.
        tx, ty = math.cos(theta) * direction, math.sin(theta) * direction
        nx, ny = -ty, tx  # normal
        tip = self._point(cx, cy, r, pos)
        length, half = 9.0, self._LANE_WIDTH * 0.7
        base_x, base_y = tip.x() - length * tx, tip.y() - length * ty
        poly = QPolygonF([
            tip,
            QPointF(base_x + half * nx, base_y + half * ny),
            QPointF(base_x - half * nx, base_y - half * ny),
        ])
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QBrush(color))
        painter.drawPolygon(poly)

    @staticmethod
    def _fmt_bp(pos: int) -> str:
        if pos == 0:
            return "1"
        if pos >= 1000 and pos % 1000 == 0:
            return f"{pos // 1000}k"
        return str(pos)
