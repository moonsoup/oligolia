"""Bundled application fonts (issue #46).

Ships real, OFL-licensed font files rather than relying on system font names
that may not exist on the user's machine (which caused Qt's slow
"Populating font family aliases" fallback and inconsistent rendering):

- Inter (UI chrome, headings via its SemiBold weight) — SIL OFL.
- JetBrains Mono (sequence / code / table monospace text) — SIL OFL.

License texts live alongside the .ttf files in ``assets/fonts/``.
"""

from __future__ import annotations

import os
import sys

FONT_UI = "Inter"
FONT_MONO = "JetBrains Mono"

_FONT_FILES = ("Inter-Regular.ttf", "Inter-SemiBold.ttf", "JetBrainsMono-Regular.ttf")


def _fonts_dir() -> str:
    # Frozen (PyInstaller) bundles assets under sys._MEIPASS; in dev the fonts
    # live at <project root>/assets/fonts (this file is gui/fonts.py).
    base = getattr(sys, "_MEIPASS", None)
    if base is None:
        base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base, "assets", "fonts")


def load_bundled_fonts() -> None:
    """Register the bundled fonts with Qt. Call once after QApplication exists."""
    from PyQt6.QtGui import QFontDatabase

    fonts_dir = _fonts_dir()
    for filename in _FONT_FILES:
        path = os.path.join(fonts_dir, filename)
        if os.path.exists(path):
            QFontDatabase.addApplicationFont(path)
