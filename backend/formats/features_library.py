"""Loader for the bundled reference library of common vector parts (issue #42).

Reads ``backend/data/common_features.json`` into structured
:class:`ReferenceFeature` records for downstream consumers (the auto-annotation
homology scan, #43).
"""

from __future__ import annotations

import json
from pathlib import Path

from ..models.sequence import ReferenceFeature

# backend/formats/features_library.py -> backend/data/common_features.json
_DATA_PATH = Path(__file__).resolve().parents[1] / "data" / "common_features.json"


def load_common_features(path: str | Path | None = None) -> list[ReferenceFeature]:
    """Return the bundled reference parts as structured records.

    A fresh list of fresh models is returned on each call so callers can freely
    mutate results without affecting the source data.
    """
    raw = json.loads(Path(path or _DATA_PATH).read_text(encoding="utf-8"))
    return [ReferenceFeature(**p) for p in raw["parts"]]
