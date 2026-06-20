"""Update available dialog — shown when a newer release is detected."""

from __future__ import annotations
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTextEdit, QProgressBar, QApplication,
)
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtCore import QUrl

from .updater import UpdateInfo, DownloadWorker, launch_installer, _platform_asset_name

try:
    from version import VERSION
except ImportError:
    VERSION = "0.0.0"


class UpdateDialog(QDialog):
    def __init__(self, info: UpdateInfo, parent=None) -> None:
        super().__init__(parent)
        self._info = info
        self._worker: DownloadWorker | None = None
        self.setWindowTitle("Update Available")
        self.setMinimumWidth(520)
        self.setMinimumHeight(360)
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setSpacing(12)
        layout.setContentsMargins(20, 20, 20, 20)

        # Header
        title = QLabel(f"Oligolia {self._info.version} is available")
        title.setObjectName("heading")
        font = title.font()
        font.setPointSize(16)
        font.setBold(True)
        title.setFont(font)
        layout.addWidget(title)

        current = QLabel(f"You have version {VERSION}.")
        current.setObjectName("subheading")
        layout.addWidget(current)

        # Release notes
        notes_label = QLabel("What's new:")
        notes_label.setObjectName("subheading")
        layout.addWidget(notes_label)

        self._notes = QTextEdit()
        self._notes.setReadOnly(True)
        self._notes.setMarkdown(self._info.body or "_No release notes provided._")
        self._notes.setMaximumHeight(180)
        layout.addWidget(self._notes)

        # Progress bar (hidden until download starts)
        self._progress = QProgressBar()
        self._progress.setRange(0, 100)
        self._progress.hide()
        layout.addWidget(self._progress)

        self._status = QLabel("")
        self._status.setObjectName("subheading")
        self._status.hide()
        layout.addWidget(self._status)

        # Buttons
        btn_row = QHBoxLayout()

        btn_skip = QPushButton("Skip This Version")
        btn_skip.clicked.connect(self.reject)
        btn_row.addWidget(btn_skip)

        btn_browser = QPushButton("Open in Browser")
        btn_browser.clicked.connect(lambda: QDesktopServices.openUrl(QUrl(self._info.html_url)))
        btn_row.addWidget(btn_browser)

        btn_row.addStretch()

        self._btn_install = QPushButton("Download & Install")
        self._btn_install.setObjectName("primary")
        self._btn_install.setMinimumWidth(160)
        self._btn_install.clicked.connect(self._start_download)
        btn_row.addWidget(self._btn_install)

        layout.addLayout(btn_row)

    def _start_download(self) -> None:
        url = self._info.download_url
        filename = _platform_asset_name(self._info.version)

        # If URL is just the releases page, open browser instead
        if not (url.endswith(".dmg") or url.endswith(".exe") or url.endswith(".AppImage")):
            QDesktopServices.openUrl(QUrl(url))
            self.accept()
            return

        self._btn_install.setEnabled(False)
        self._btn_install.setText("Downloading…")
        self._progress.show()
        self._progress.setValue(0)
        self._status.show()
        self._status.setText(f"Downloading {filename}…")

        self._worker = DownloadWorker(url, filename)
        self._worker.progress.connect(self._progress.setValue)
        self._worker.finished.connect(self._on_downloaded)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_downloaded(self, path: str) -> None:
        self._progress.setValue(100)
        self._status.setText("Launching installer…")
        launch_installer(path)
        # Quit the app so the installer can replace it
        QApplication.quit()

    def _on_error(self, err: str) -> None:
        self._status.setText(f"Download failed: {err}")
        self._btn_install.setEnabled(True)
        self._btn_install.setText("Download & Install")
        self._progress.hide()
