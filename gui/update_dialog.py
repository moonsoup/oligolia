"""Update available dialog — patch or full installer, with progress."""

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

from .updater import (
    UpdateInfo, DownloadWorker,
    apply_code_patch, launch_full_installer, restart_app,
)

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
        self.setMinimumWidth(540)
        self.setMinimumHeight(380)
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setSpacing(12)
        layout.setContentsMargins(20, 20, 20, 20)

        # Header
        title = QLabel(f"Oligolia {self._info.version} is available")
        font = title.font()
        font.setPointSize(16)
        font.setBold(True)
        title.setFont(font)
        title.setStyleSheet("color: #4ade80;")
        layout.addWidget(title)

        layout.addWidget(QLabel(f"You have version {VERSION}."))

        # Update type badge
        if self._info.can_patch:
            badge = QLabel("⚡ Code patch — only downloads what changed (~11 MB)")
            badge.setStyleSheet("color: #4ade80; font-size: 12px;")
        else:
            badge = QLabel("📦 Full installer (~260 MB) — runtime dependencies changed")
            badge.setStyleSheet("color: #fbbf24; font-size: 12px;")
        layout.addWidget(badge)

        # Release notes
        layout.addWidget(QLabel("What's new:"))
        self._notes = QTextEdit()
        self._notes.setReadOnly(True)
        self._notes.setMarkdown(self._info.body or "_No release notes provided._")
        self._notes.setMaximumHeight(170)
        layout.addWidget(self._notes)

        # Progress
        self._progress = QProgressBar()
        self._progress.setRange(0, 100)
        self._progress.hide()
        layout.addWidget(self._progress)

        self._status = QLabel("")
        self._status.setStyleSheet("color: #94a3b8; font-size: 12px;")
        self._status.hide()
        layout.addWidget(self._status)

        # Buttons
        btn_row = QHBoxLayout()
        btn_skip = QPushButton("Skip")
        btn_skip.clicked.connect(self.reject)
        btn_row.addWidget(btn_skip)

        btn_browser = QPushButton("Open in Browser")
        btn_browser.clicked.connect(
            lambda: QDesktopServices.openUrl(QUrl(self._info.html_url))
        )
        btn_row.addWidget(btn_browser)
        btn_row.addStretch()

        label = "Download Patch (~11 MB)" if self._info.can_patch else "Download & Install (~260 MB)"
        self._btn_install = QPushButton(label)
        self._btn_install.setObjectName("primary")
        self._btn_install.setMinimumWidth(200)
        self._btn_install.clicked.connect(self._start_download)
        btn_row.addWidget(self._btn_install)
        layout.addLayout(btn_row)

    def _start_download(self) -> None:
        url = self._info.download_url
        filename = url.split("/")[-1] if "/" in url else "update"

        # Fallback: open browser if no direct download URL
        if not any(url.endswith(ext) for ext in (".dmg", ".exe", ".AppImage", ".tar.gz", ".zip")):
            QDesktopServices.openUrl(QUrl(url))
            self.accept()
            return

        self._btn_install.setEnabled(False)
        self._btn_install.setText("Downloading…")
        self._progress.show()
        self._progress.setValue(0)
        self._status.show()
        self._status.setText(f"Downloading {filename} ({self._info.download_size_hint})…")

        self._worker = DownloadWorker(url, filename)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_downloaded)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_progress(self, pct: int) -> None:
        if pct == -1:
            # Server sent no Content-Length — switch to indeterminate (pulsing) mode
            self._progress.setRange(0, 0)
        else:
            if self._progress.maximum() == 0:
                self._progress.setRange(0, 100)
            self._progress.setValue(pct)

    def _on_downloaded(self, path: str) -> None:
        self._progress.setValue(100)
        if self._info.can_patch:
            self._status.setText("Applying patch…")
            QApplication.processEvents()
            try:
                apply_code_patch(path)
                self._status.setText("Patch applied — restarting…")
                QApplication.processEvents()
                QApplication.quit()
                restart_app()
            except Exception as e:
                self._status.setText(f"Patch failed: {e} — falling back to full installer")
                self._btn_install.setText("Download Full Installer")
                self._btn_install.setEnabled(True)
                # Update info to force full install next click
                self._info.requires_full = True
        else:
            self._status.setText("Launching installer…")
            launch_full_installer(path)
            QApplication.quit()

    def _on_error(self, err: str) -> None:
        self._status.setText(f"Download failed: {err}")
        self._btn_install.setEnabled(True)
        self._btn_install.setText("Retry")
        self._progress.hide()
