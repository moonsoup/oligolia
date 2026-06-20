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
    app.setApplicationName("Oligolia")
    app.setApplicationVersion("0.1.0")
    app.setOrganizationName("Oligolia")

    # Use a clean system font
    font = QFont()
    font.setFamily("-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif")
    font.setPointSize(13)
    app.setFont(font)

    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
