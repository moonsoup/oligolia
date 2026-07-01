"""Shared feature-type → color mapping.

Single source of truth for how a feature type (CDS, exon, promoter, …) is
colored, so the sequence highlighter, the feature table, the minimap, and the
plasmid map all agree — a given feature type looks the same everywhere.
"""

from __future__ import annotations

from PyQt6.QtGui import QColor

from backend.models.sequence import Annotation

# Distinct dark-theme background hues, cycled across feature types.
FEATURE_PALETTE: list[str] = [
    "#1e3a5f", "#5f1e3a", "#3a5f1e", "#5f4a1e",
    "#1e5f4a", "#4a1e5f", "#5f1e1e", "#274060",
]


def feature_color_map(annotations: list[Annotation]) -> dict[str, QColor]:
    """Map each feature type present to a stable :class:`QColor`.

    Types are sorted by name so the assignment is deterministic regardless of
    annotation order — the same set of feature types always yields the same
    colors across every view.
    """
    types = sorted({a.feature_type for a in (annotations or [])})
    return {
        t: QColor(FEATURE_PALETTE[i % len(FEATURE_PALETTE)])
        for i, t in enumerate(types)
    }
