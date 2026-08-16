"""Detects the optional, separately-installed "Oligolia Structure Viewer" companion
app (structure_viewer/) — a small PyQt6-WebEngine app kept out of the core Oligolia
installer to avoid bundling a Chromium runtime into every user's download.

Detection order: (1) OS-standard install location, (2) a user-provided path cached
in ~/.oligolia/config.json (set once via a "Locate Structure Viewer…" file picker
in the calling panel, if auto-detection fails). No silent installs, no auto-download
— the user installs the companion app themselves.
"""

from __future__ import annotations
import json
import platform
from pathlib import Path

CONFIG_PATH = Path.home() / ".oligolia" / "config.json"

_STANDARD_LOCATIONS = {
    "Darwin": [Path("/Applications/Oligolia Structure Viewer.app/Contents/MacOS/Oligolia Structure Viewer")],
    "Windows": [Path(r"C:\Program Files\Oligolia Structure Viewer\Oligolia Structure Viewer.exe")],
    "Linux": [
        Path.home() / ".local" / "share" / "oligolia" / "structure-viewer" / "Oligolia-Structure-Viewer.AppImage",
        Path("/opt/oligolia-structure-viewer/Oligolia-Structure-Viewer.AppImage"),
    ],
}


def _cached_path() -> Path | None:
    if not CONFIG_PATH.exists():
        return None
    try:
        data = json.loads(CONFIG_PATH.read_text())
    except (json.JSONDecodeError, OSError):
        return None
    path = data.get("structure_viewer_path")
    return Path(path) if path else None


def set_cached_path(path: Path) -> None:
    """Persist a user-provided viewer location (e.g. from a file picker)."""
    CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    data = {}
    if CONFIG_PATH.exists():
        try:
            data = json.loads(CONFIG_PATH.read_text())
        except (json.JSONDecodeError, OSError):
            data = {}
    data["structure_viewer_path"] = str(path)
    CONFIG_PATH.write_text(json.dumps(data))


def find_structure_viewer() -> Path | None:
    """Return the companion viewer's executable path, or None if not installed."""
    for candidate in _STANDARD_LOCATIONS.get(platform.system(), []):
        if candidate.exists():
            return candidate

    cached = _cached_path()
    if cached is not None and cached.exists():
        return cached

    return None
