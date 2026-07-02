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

from PyQt6.QtCore import QPointF, QRectF, Qt, pyqtSignal
from PyQt6.QtGui import QBrush, QColor, QFont, QPainter, QPen, QPolygonF
from PyQt6.QtWidgets import QWidget

from backend.models.sequence import Annotation, Sequence, Strand
from backend.routers.primers import RestrictionRequest, restriction_sites
from gui.panels.feature_colors import feature_color_map

_CUT_COLOR = "#b91c1c"  # distinct red for enzyme cut ticks/labels

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

    feature_clicked = pyqtSignal(int, int)  # (start, end) of a clicked arc

    _LANE_WIDTH = 9  # arc thickness / lane spacing (px)

    def __init__(self) -> None:
        super().__init__()
        self._name = ""
        self._length = 0
        self._is_circular = False
        self._lanes: list[tuple[Annotation, int]] = []
        self._color_map: dict[str, QColor] = {}
        self._sites: list = []  # RestrictionSite results
        self._unique_only = True  # SnapGene default: show only single cutters
        self._highlight: tuple[int, int] | None = None  # selected span to glow
        self.setMinimumSize(240, 240)

    def set_highlight(self, span: tuple[int, int] | None) -> None:
        """Glow the arc(s) overlapping [start, end); None clears it."""
        self._highlight = span
        self.update()

    def _ring_radius(self) -> float:
        return min(self.width(), self.height()) / 2 - 50  # margin for labels

    def _lane_radius(self, lane: int, ring: float) -> float:
        return ring - 12 - lane * (self._LANE_WIDTH + 2)

    def _overlaps_highlight(self, ann: Annotation) -> bool:
        if self._highlight is None:
            return False
        hs, he = self._highlight
        return not (ann.end <= hs or ann.start >= he)

    def mousePressEvent(self, event) -> None:  # noqa: N802 (Qt override)
        """Hit-test feature arcs; emit the clicked feature's span (issue #26)."""
        if not self._lanes or not self._length:
            return
        ring = self._ring_radius()
        dx = event.position().x() - self.width() / 2
        dy = event.position().y() - self.height() / 2
        r_click = math.hypot(dx, dy)
        theta = math.atan2(dx, -dy)  # 0 at top (12 o'clock), clockwise positive
        if theta < 0:
            theta += 2 * math.pi
        pos = theta / (2 * math.pi) * self._length

        best, best_dr = None, 1e9
        for ann, lane in self._lanes:
            r_lane = self._lane_radius(lane, ring)
            dr = abs(r_click - r_lane)
            if dr <= self._LANE_WIDTH and ann.start <= pos <= ann.end and dr < best_dr:
                best, best_dr = ann, dr
        if best is not None:
            self.feature_clicked.emit(best.start, best.end)

    def set_sequence(self, seq: Sequence) -> None:
        """Point the map at a sequence (intended for circular ones)."""
        self._name = seq.name or seq.id
        self._length = seq.length
        self._is_circular = seq.is_circular
        self._lanes = _assign_lanes(seq.annotations)
        self._color_map = feature_color_map(seq.annotations)
        # Reuse the backend restriction search (circular-aware, #21) in-process.
        try:
            self._sites = restriction_sites(
                RestrictionRequest(template=seq.seq, is_circular=seq.is_circular))
        except Exception:
            self._sites = []
        self.update()

    def set_unique_cutters_only(self, unique_only: bool) -> None:
        """Toggle between unique cutters (default) and all restriction sites."""
        self._unique_only = unique_only
        self.update()

    def _point(self, cx: float, cy: float, radius: float, pos: int) -> QPointF:
        theta = 2 * math.pi * (pos / self._length) if self._length else 0.0
        return QPointF(cx + radius * math.sin(theta), cy - radius * math.cos(theta))

    def paintEvent(self, event) -> None:  # noqa: N802 (Qt override)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        w, h = self.width(), self.height()
        cx, cy = w / 2, h / 2
        radius = self._ring_radius()  # margin for bp-ruler + enzyme labels

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
        label_font = QFont("JetBrains Mono", 8)
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

        # ── Restriction cut sites (issue #25) ─────────────────────────────
        self._draw_cut_sites(painter, cx, cy, radius)

        # ── Centre label: name + size + topology ──────────────────────────
        painter.setPen(QPen(QColor("#e2e8f0"), 1))
        painter.setFont(QFont("Inter", 11, QFont.Weight.Bold))
        name_rect = QRectF(cx - radius, cy - 16, 2 * radius, 18)
        painter.drawText(name_rect, Qt.AlignmentFlag.AlignCenter, self._name or "(unnamed)")
        painter.setPen(QPen(QColor("#94a3b8"), 1))
        painter.setFont(QFont("Inter", 9))
        topo = "circular" if self._is_circular else "linear"
        size_rect = QRectF(cx - radius, cy + 4, 2 * radius, 16)
        painter.drawText(size_rect, Qt.AlignmentFlag.AlignCenter, f"{self._length:,} bp · {topo}")
        painter.end()

    def _draw_features(self, painter: QPainter, cx: float, cy: float, radius: float) -> None:
        """Draw each annotation as a colored arc + strand arrowhead + label."""
        label_font = QFont("Inter", 8)
        for ann, lane in self._lanes:
            r = self._lane_radius(lane, radius)
            if r < radius * 0.35:  # too many overlapping lanes to fit; skip inner ones
                continue
            color = self._color_map.get(ann.feature_type, QColor("#94a3b8"))

            # Arc: Qt angles are 1/16 deg, 0 at 3 o'clock, CCW positive. Our bp 0
            # is at 12 o'clock increasing clockwise, so map bp p -> 90 - 360*p/L.
            span_bp = max(ann.end - ann.start, 1)
            start_angle = (90.0 - 360.0 * ann.start / self._length) * 16
            span_angle = -(360.0 * span_bp / self._length) * 16
            rect = QRectF(cx - r, cy - r, 2 * r, 2 * r)
            # Selection glow (issue #26): a wide translucent halo under arcs that
            # overlap the current sequence-view selection.
            if self._overlaps_highlight(ann):
                painter.setPen(QPen(QColor(250, 204, 21, 160), self._LANE_WIDTH + 8,
                                    Qt.PenStyle.SolidLine, Qt.PenCapStyle.FlatCap))
                painter.drawArc(rect, round(start_angle), round(span_angle))
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

    def _draw_cut_sites(self, painter: QPainter, cx: float, cy: float, radius: float) -> None:
        """Draw restriction cut-site ticks + enzyme labels on the outer edge.

        Enzyme labels are stacked in left/right side columns with leader lines
        back to their ticks, so densely-clustered cutters (e.g. a plasmid MCS)
        stay legible instead of drawing over each other.
        """
        display: list[tuple[str, int]] = []
        for site in self._sites:
            if self._unique_only and site.count != 1:
                continue
            for pos in site.positions:
                display.append((site.enzyme, pos))
        if not display:
            return

        color = QColor(_CUT_COLOR)
        painter.setFont(QFont("Inter", 8))

        # 1) Radial ticks at each cut position.
        painter.setPen(QPen(color, 1))
        right: list[tuple[float, str, QPointF]] = []
        left: list[tuple[float, str, QPointF]] = []
        for enzyme, pos in display:
            painter.drawLine(self._point(cx, cy, radius, pos),
                             self._point(cx, cy, radius + 9, pos))
            tip = self._point(cx, cy, radius + 9, pos)
            (right if tip.x() >= cx else left).append((tip.y(), enzyme, tip))

        # 2) Side-stacked labels with leader lines (vertical de-collision).
        gap = 14
        for side, column in (("R", right), ("L", left)):
            column.sort(key=lambda t: t[0])
            col_x = cx + radius + 18 if side == "R" else cx - radius - 18
            last_y = -1e9
            for _, enzyme, tip in column:
                y = max(tip.y(), last_y + gap)
                last_y = y
                painter.setPen(QPen(color, 1))
                painter.drawLine(tip, QPointF(col_x, y))
                if side == "R":
                    box = QRectF(col_x + 2, y - 8, 120, 16)
                    align = Qt.AlignmentFlag.AlignLeft
                else:
                    box = QRectF(col_x - 122, y - 8, 120, 16)
                    align = Qt.AlignmentFlag.AlignRight
                painter.drawText(box, align | Qt.AlignmentFlag.AlignVCenter, enzyme)

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
