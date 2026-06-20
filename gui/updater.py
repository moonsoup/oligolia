"""
Auto-update checker and installer launcher for Oligolia.

Flow:
  1. On startup: check GitHub Releases API in a background thread
  2. If a newer version exists: show UpdateDialog (non-blocking)
  3. User clicks "Download & Install": downloads platform installer to temp dir
  4. Launches the installer, then quits (installer replaces the app)

Platform installer files expected on GitHub Releases:
  macOS  : Oligolia-{version}-mac.dmg
  Windows: Oligolia-{version}-Setup.exe
  Linux  : Oligolia-{version}-x86_64.AppImage
"""

from __future__ import annotations
import sys
import os
import platform
import tempfile
import subprocess
from packaging.version import Version

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import httpx
from PyQt6.QtCore import QThread, pyqtSignal

try:
    from version import VERSION, APP_NAME, RELEASES_API_URL, RELEASES_PAGE_URL
except ImportError:
    VERSION = "0.0.0"
    APP_NAME = "Oligolia"
    RELEASES_API_URL = "https://api.github.com/repos/moonsoup/oligolia/releases/latest"
    RELEASES_PAGE_URL = "https://github.com/moonsoup/oligolia/releases"


def _platform_asset_name(version: str) -> str:
    """Return the expected release asset filename for this platform."""
    system = platform.system()
    arch = platform.machine().lower()
    if system == "Darwin":
        return f"Oligolia-{version}-mac.dmg"
    if system == "Windows":
        return f"Oligolia-{version}-Setup.exe"
    # Linux
    arch_tag = "x86_64" if arch in ("x86_64", "amd64") else arch
    return f"Oligolia-{version}-{arch_tag}.AppImage"


class UpdateInfo:
    def __init__(self, version: str, body: str, download_url: str, html_url: str) -> None:
        self.version = version
        self.body = body
        self.download_url = download_url
        self.html_url = html_url

    @property
    def is_newer(self) -> bool:
        try:
            return Version(self.version) > Version(VERSION)
        except Exception:
            return False


class UpdateChecker(QThread):
    """Background thread that checks for a new release."""
    update_available = pyqtSignal(object)   # emits UpdateInfo
    check_failed = pyqtSignal(str)          # emits error message (silently ignored in UI)

    def run(self) -> None:
        try:
            with httpx.Client(timeout=10, follow_redirects=True) as client:
                r = client.get(
                    RELEASES_API_URL,
                    headers={"Accept": "application/vnd.github+json", "X-GitHub-Api-Version": "2022-11-28"},
                )
                if r.status_code == 404:
                    # Private repo or no releases yet — silently skip
                    return
                r.raise_for_status()
                data = r.json()

            tag = data.get("tag_name", "").lstrip("v")
            body = data.get("body", "")
            html_url = data.get("html_url", RELEASES_PAGE_URL)
            assets = data.get("assets", [])

            wanted_name = _platform_asset_name(tag)
            download_url = ""
            for asset in assets:
                if asset.get("name") == wanted_name:
                    download_url = asset.get("browser_download_url", "")
                    break

            if not download_url:
                download_url = html_url  # fallback: open releases page

            info = UpdateInfo(version=tag, body=body, download_url=download_url, html_url=html_url)
            if info.is_newer:
                self.update_available.emit(info)

        except Exception as e:
            self.check_failed.emit(str(e))


class DownloadWorker(QThread):
    """Downloads the installer file in the background."""
    progress = pyqtSignal(int)    # 0–100
    finished = pyqtSignal(str)    # path to downloaded file
    error = pyqtSignal(str)

    def __init__(self, url: str, filename: str) -> None:
        super().__init__()
        self._url = url
        self._filename = filename

    def run(self) -> None:
        try:
            dest = os.path.join(tempfile.gettempdir(), self._filename)
            with httpx.Client(timeout=300, follow_redirects=True) as client:
                with client.stream("GET", self._url) as r:
                    r.raise_for_status()
                    total = int(r.headers.get("content-length", 0))
                    downloaded = 0
                    with open(dest, "wb") as f:
                        for chunk in r.iter_bytes(chunk_size=65536):
                            f.write(chunk)
                            downloaded += len(chunk)
                            if total:
                                self.progress.emit(int(downloaded / total * 100))
            self.finished.emit(dest)
        except Exception as e:
            self.error.emit(str(e))


def launch_installer(path: str) -> None:
    """Open the downloaded installer using the platform-native method."""
    system = platform.system()
    if system == "Darwin":
        subprocess.Popen(["open", path])
    elif system == "Windows":
        os.startfile(path)  # type: ignore[attr-defined]
    else:
        # Linux: make executable and run
        os.chmod(path, 0o755)
        subprocess.Popen([path])
