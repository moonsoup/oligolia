"""
Auto-update checker and installer/patcher for Oligolia.

Two-tier update strategy:
  TIER 1 — Code patch (~11 MB, 96% smaller):
    Downloads only the changed executable + version file.
    Applies in-place to the installed .app / install dir.
    Used for all normal releases (bug fixes, new features).

  TIER 2 — Full installer (~260 MB):
    Full DMG / Setup.exe / AppImage.
    Used when Python runtime, PyQt6, or Biopython version changes.
    Indicated by requires_full=true in the release manifest.

Each GitHub Release includes:
  - Oligolia-{ver}-manifest.json          (always, tiny)
  - Oligolia-{ver}-mac-patch.tar.gz       (~11 MB)
  - Oligolia-{ver}-mac.dmg               (~260 MB, full)
  - Oligolia-{ver}-Setup.exe             (Windows full)
  - Oligolia-{ver}-x86_64.AppImage       (Linux full)
"""

from __future__ import annotations
import sys
import os
import platform
import tempfile
import subprocess
import tarfile
import shutil
import logging
from pathlib import Path
from packaging.version import Version

# Write update diagnostics to a log file the user can inspect
_log_path = Path(tempfile.gettempdir()) / "oligolia_update.log"
logging.basicConfig(
    filename=str(_log_path),
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)s %(message)s",
)
_log = logging.getLogger("oligolia.updater")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import httpx
from PyQt6.QtCore import QThread, pyqtSignal

try:
    from version import VERSION, RELEASES_API_URL, RELEASES_PAGE_URL
except ImportError:
    VERSION = "0.0.0"
    RELEASES_API_URL = "https://api.github.com/repos/moonsoup/oligolia/releases/latest"
    RELEASES_PAGE_URL = "https://github.com/moonsoup/oligolia/releases"


def _app_bundle_path() -> Path | None:
    """Return the path to the running .app bundle (macOS only)."""
    exe = Path(sys.executable)
    # PyInstaller: sys.executable = .app/Contents/MacOS/Oligolia
    for parent in exe.parents:
        if parent.suffix == ".app":
            return parent
    return None


def _install_dir() -> Path | None:
    """Return the directory containing the running executable (Windows/Linux)."""
    return Path(sys.executable).parent if getattr(sys, "frozen", False) else None


def _platform_patch_asset(version: str) -> str:
    if platform.system() == "Darwin":
        return f"Oligolia-{version}-mac-patch.tar.gz"
    return ""  # patch not yet implemented for Win/Linux (use full installer)


def _platform_full_asset(version: str) -> str:
    system = platform.system()
    arch = platform.machine().lower()
    if system == "Darwin":
        return f"Oligolia-{version}-mac.dmg"
    if system == "Windows":
        return f"Oligolia-{version}-Setup.exe"
    arch_tag = "x86_64" if arch in ("x86_64", "amd64") else arch
    return f"Oligolia-{version}-{arch_tag}.AppImage"


class UpdateInfo:
    def __init__(self, version: str, body: str, html_url: str,
                 patch_url: str, full_url: str, requires_full: bool,
                 min_compatible_base: str) -> None:
        self.version = version
        self.body = body
        self.html_url = html_url
        self.patch_url = patch_url          # empty string = not available
        self.full_url = full_url
        self.requires_full = requires_full
        self.min_compatible_base = min_compatible_base

    @property
    def is_newer(self) -> bool:
        try:
            return Version(self.version) > Version(VERSION)
        except Exception:
            return False

    @property
    def can_patch(self) -> bool:
        """True if a code patch is available and compatible with this installation."""
        if self.requires_full or not self.patch_url:
            return False
        try:
            return Version(VERSION) >= Version(self.min_compatible_base)
        except Exception:
            return False

    @property
    def download_url(self) -> str:
        return self.patch_url if self.can_patch else self.full_url

    @property
    def download_size_hint(self) -> str:
        return "~11 MB (code patch)" if self.can_patch else "~260 MB (full installer)"


class UpdateChecker(QThread):
    """Background thread — checks GitHub Releases API for a newer version."""
    update_available = pyqtSignal(object)
    check_failed = pyqtSignal(str)

    def run(self) -> None:
        _log.info("Update check started. Current VERSION=%s", VERSION)
        try:
            with httpx.Client(timeout=10, follow_redirects=True) as client:
                r = client.get(
                    RELEASES_API_URL,
                    headers={"Accept": "application/vnd.github+json",
                             "X-GitHub-Api-Version": "2022-11-28"},
                )
                _log.info("GitHub API status: %s", r.status_code)
                if r.status_code == 404:
                    _log.info("No releases found (404) — skipping")
                    return
                r.raise_for_status()
                data = r.json()

            tag = data.get("tag_name", "").lstrip("v")
            body = data.get("body", "")
            html_url = data.get("html_url", RELEASES_PAGE_URL)
            assets: dict[str, str] = {
                a["name"]: a["browser_download_url"]
                for a in data.get("assets", [])
            }
            _log.info("Latest tag: %s  Assets: %s", tag, list(assets.keys()))

            # Try to fetch the manifest for patch/full decision
            manifest_name = f"Oligolia-{tag}-manifest.json"
            manifest_url = assets.get(manifest_name, "")
            requires_full = False
            min_compat = "0.0.0"
            patch_name = _platform_patch_asset(tag)
            full_name = _platform_full_asset(tag)
            _log.info("Looking for manifest: %s  found=%s", manifest_name, bool(manifest_url))

            if manifest_url:
                try:
                    with httpx.Client(timeout=10, follow_redirects=True) as mclient:
                        mresp = mclient.get(manifest_url)
                    manifest = mresp.json()
                    requires_full = manifest.get("requires_full", False)
                    min_compat = manifest.get("min_compatible_base", "0.0.0")
                    m_assets = manifest.get("assets", {})
                    if "darwin_patch" in m_assets and platform.system() == "Darwin":
                        patch_name = m_assets["darwin_patch"]
                    if "darwin_full" in m_assets and platform.system() == "Darwin":
                        full_name = m_assets["darwin_full"]
                    _log.info("Manifest: requires_full=%s min_compat=%s patch_name=%s",
                              requires_full, min_compat, patch_name)
                except Exception as e:
                    _log.warning("Manifest fetch failed: %s", e)

            patch_url = assets.get(patch_name, "")
            full_url = assets.get(full_name, html_url)
            _log.info("patch_url=%s", patch_url or "(none)")
            _log.info("full_url=%s", full_url)

            info = UpdateInfo(
                version=tag, body=body, html_url=html_url,
                patch_url=patch_url, full_url=full_url,
                requires_full=requires_full, min_compatible_base=min_compat,
            )
            _log.info("is_newer=%s can_patch=%s download_url=%s",
                      info.is_newer, info.can_patch, info.download_url)
            if info.is_newer:
                self.update_available.emit(info)

        except Exception as e:
            _log.error("Update check failed: %s", e)
            self.check_failed.emit(str(e))


class DownloadWorker(QThread):
    """Downloads an update file with streaming progress."""
    progress = pyqtSignal(int)   # 0–100
    finished = pyqtSignal(str)   # local file path
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
                    if not total:
                        # Signal indeterminate mode — dialog will pulse the bar
                        self.progress.emit(-1)
                    with open(dest, "wb") as f:
                        for chunk in r.iter_bytes(chunk_size=65536):
                            f.write(chunk)
                            downloaded += len(chunk)
                            if total:
                                self.progress.emit(int(downloaded / total * 100))
            self.finished.emit(dest)
        except Exception as e:
            self.error.emit(str(e))


def apply_code_patch(patch_path: str) -> None:
    """
    Apply a macOS code patch (tar.gz) to the running .app bundle.

    The patch contains:
      Contents/MacOS/Oligolia        ← replace the main executable
      Contents/Resources/version.py  ← update version number
      Contents/Info.plist            ← update plist version strings
    """
    app = _app_bundle_path()
    if not app:
        raise RuntimeError("Cannot locate .app bundle — are you running from the installed app?")

    with tarfile.open(patch_path, "r:gz") as tar:
        for member in tar.getmembers():
            # Resolve destination and reject any path that escapes the .app bundle
            dest = (app / member.name).resolve()
            if not str(dest).startswith(str(app.resolve())):
                raise RuntimeError(
                    f"Patch contains unsafe path '{member.name}' — aborting"
                )
            dest.parent.mkdir(parents=True, exist_ok=True)
            if member.isfile():
                with tar.extractfile(member) as src, open(dest, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                # Preserve executable bit for the main binary
                if "MacOS/" in member.name:
                    os.chmod(dest, 0o755)

    # Clear code signature — patched binary won't match original sig
    subprocess.run(["codesign", "--remove-signature", str(app)],
                   capture_output=True)

    # Re-sign ad-hoc so Gatekeeper accepts it
    subprocess.run(["codesign", "--force", "--deep", "--sign", "-", str(app)],
                   capture_output=True)


def launch_full_installer(path: str) -> None:
    """Launch a full platform installer."""
    system = platform.system()
    if system == "Darwin":
        subprocess.Popen(["open", path])
    elif system == "Windows":
        os.startfile(path)  # type: ignore[attr-defined]
    else:
        os.chmod(path, 0o755)
        subprocess.Popen([path])


def restart_app() -> None:
    """Restart the current process (after a patch has been applied)."""
    exe = sys.executable
    os.execv(exe, [exe] + sys.argv)
