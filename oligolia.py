#!/usr/bin/env python3
"""
Oligolia — Standalone gene editing platform.

Usage:
    python oligolia.py              # Launch the GUI
    python oligolia.py --help       # Show options
"""

import sys
import os

# Ensure the project root is on the path so both `backend` and `gui` import
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QFont

from gui.main_window import MainWindow


def main() -> None:
    app = QApplication(sys.argv)
    from version import VERSION
    app.setApplicationName("Oligolia")
    app.setApplicationVersion(VERSION)
    app.setOrganizationName("Oligolia")

    # Register bundled fonts so rendering is consistent regardless of what is
    # installed on the user's machine (issue #46), then use Inter for UI chrome.
    from gui.fonts import FONT_UI, load_bundled_fonts
    load_bundled_fonts()
    app.setFont(QFont(FONT_UI, 13))

    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
