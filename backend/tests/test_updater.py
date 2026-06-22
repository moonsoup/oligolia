"""
Tests for the two-tier update system (gui/updater.py).
Uses unittest.mock to simulate GitHub API responses and filesystem state
so these tests run without a network connection or a PyInstaller bundle.
"""

import io
import os
import sys
import tarfile
import tempfile
import platform
from pathlib import Path
from unittest.mock import MagicMock, patch, PropertyMock

import pytest

# Updater lives in gui/, not backend/ — add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from gui.updater import (
    UpdateInfo,
    _platform_patch_asset,
    _platform_full_asset,
    _app_bundle_path,
    apply_code_patch,
    restart_app,
)

CURRENT = "0.3.0"


# ── UpdateInfo.is_newer ───────────────────────────────────────────────────────

def make_info(version, requires_full=False, patch_url="", full_url="http://x/full.dmg",
              min_compat="0.0.0"):
    return UpdateInfo(
        version=version, body="", html_url="http://x",
        patch_url=patch_url, full_url=full_url,
        requires_full=requires_full, min_compatible_base=min_compat,
    )


def test_is_newer_true():
    assert make_info("0.4.0").is_newer is True

def test_is_newer_false_same():
    assert make_info(CURRENT).is_newer is False

def test_is_newer_false_older():
    assert make_info("0.2.9").is_newer is False

def test_is_newer_with_v_prefix():
    # tag_name already has "v" stripped in UpdateChecker, but guard against it
    info = make_info("v0.4.0")
    # packaging.version.Version("v0.4.0") raises InvalidVersion
    # is_newer must return False (not crash) when version is malformed
    result = info.is_newer
    assert isinstance(result, bool)

def test_is_newer_malformed_version_does_not_crash():
    assert make_info("not-a-version").is_newer is False


# ── UpdateInfo.can_patch ─────────────────────────────────────────────────────

def test_can_patch_true():
    info = make_info("0.4.0", patch_url="http://x/patch.tar.gz", min_compat="0.3.0")
    assert info.can_patch is True

def test_can_patch_false_requires_full():
    info = make_info("0.4.0", patch_url="http://x/patch.tar.gz",
                     requires_full=True, min_compat="0.0.0")
    assert info.can_patch is False

def test_can_patch_false_no_patch_url():
    info = make_info("0.4.0", patch_url="", min_compat="0.0.0")
    assert info.can_patch is False

def test_can_patch_false_base_too_old():
    # Current version 0.3.0, patch requires base >= 0.3.5 → can't patch
    info = make_info("0.4.0", patch_url="http://x/patch.tar.gz", min_compat="0.3.5")
    assert info.can_patch is False

def test_can_patch_false_malformed_min_compat():
    info = make_info("0.4.0", patch_url="http://x/patch.tar.gz", min_compat="banana")
    assert info.can_patch is False


# ── UpdateInfo.download_url ──────────────────────────────────────────────────

def test_download_url_uses_patch_when_eligible():
    info = make_info("0.4.0", patch_url="http://x/patch.tar.gz",
                     full_url="http://x/full.dmg", min_compat="0.0.0")
    assert info.download_url == "http://x/patch.tar.gz"

def test_download_url_uses_full_when_no_patch():
    info = make_info("0.4.0", patch_url="", full_url="http://x/full.dmg")
    assert info.download_url == "http://x/full.dmg"

def test_download_url_uses_full_when_requires_full():
    info = make_info("0.4.0", patch_url="http://x/patch.tar.gz",
                     full_url="http://x/full.dmg", requires_full=True)
    assert info.download_url == "http://x/full.dmg"

def test_download_url_never_empty_when_full_url_set():
    info = make_info("0.4.0", patch_url="", full_url="http://x/release-page")
    assert info.download_url != ""


# ── Asset name helpers ───────────────────────────────────────────────────────

def test_patch_asset_name_darwin():
    with patch("platform.system", return_value="Darwin"):
        assert _platform_patch_asset("0.4.0") == "Oligolia-0.4.0-mac-patch.tar.gz"

def test_patch_asset_name_windows_empty():
    # Patch not yet implemented for Windows
    with patch("platform.system", return_value="Windows"):
        assert _platform_patch_asset("0.4.0") == ""

def test_full_asset_name_darwin():
    with patch("platform.system", return_value="Darwin"):
        assert _platform_full_asset("0.4.0") == "Oligolia-0.4.0-mac.dmg"

def test_full_asset_name_windows():
    with patch("platform.system", return_value="Windows"):
        assert _platform_full_asset("0.4.0") == "Oligolia-0.4.0-Setup.exe"

def test_full_asset_name_linux_x86():
    with patch("platform.system", return_value="Linux"), \
         patch("platform.machine", return_value="x86_64"):
        assert _platform_full_asset("0.4.0") == "Oligolia-0.4.0-x86_64.AppImage"


# ── UpdateChecker API parsing ─────────────────────────────────────────────────

MOCK_RELEASE_WITH_ASSETS = {
    "tag_name": "0.4.0",
    "body": "Bug fixes",
    "html_url": "https://github.com/moonsoup/oligolia/releases/tag/0.4.0",
    "assets": [
        {"name": "Oligolia-0.4.0-mac-patch.tar.gz",
         "browser_download_url": "https://github.com/.../mac-patch.tar.gz"},
        {"name": "Oligolia-0.4.0-mac.dmg",
         "browser_download_url": "https://github.com/.../mac.dmg"},
        {"name": "Oligolia-0.4.0-manifest.json",
         "browser_download_url": "https://github.com/.../manifest.json"},
    ],
}

MOCK_MANIFEST = {
    "requires_full": False,
    "min_compatible_base": "0.3.0",
}

MOCK_RELEASE_NO_ASSETS = {
    "tag_name": "0.4.0",
    "body": "Bug fixes",
    "html_url": "https://github.com/moonsoup/oligolia/releases/tag/0.4.0",
    "assets": [],
}


def _mock_http_response(json_body, status=200):
    resp = MagicMock()
    resp.status_code = status
    resp.json.return_value = json_body
    resp.raise_for_status = MagicMock()
    return resp


def test_update_checker_emits_when_newer():
    """UpdateChecker should emit update_available signal when a newer version exists."""
    from gui.updater import UpdateChecker

    emitted = []
    checker = UpdateChecker()
    checker.update_available = MagicMock()
    checker.update_available.emit = lambda info: emitted.append(info)
    checker.check_failed = MagicMock()

    manifest_resp = _mock_http_response(MOCK_MANIFEST)
    release_resp = _mock_http_response(MOCK_RELEASE_WITH_ASSETS)

    with patch("platform.system", return_value="Darwin"), \
         patch("httpx.Client") as mock_client_cls:
        client_instance = MagicMock()
        mock_client_cls.return_value.__enter__.return_value = client_instance
        client_instance.get.side_effect = [release_resp, manifest_resp]

        checker.run()

    assert len(emitted) == 1
    info = emitted[0]
    assert info.version == "0.4.0"
    assert info.can_patch is True
    assert info.patch_url == "https://github.com/.../mac-patch.tar.gz"


def test_update_checker_silent_when_same_version():
    """No signal emitted when latest release == current version."""
    from gui.updater import UpdateChecker

    emitted = []
    checker = UpdateChecker()
    checker.update_available = MagicMock()
    checker.update_available.emit = lambda info: emitted.append(info)
    checker.check_failed = MagicMock()

    same_release = {**MOCK_RELEASE_WITH_ASSETS, "tag_name": CURRENT}
    release_resp = _mock_http_response(same_release)
    manifest_resp = _mock_http_response(MOCK_MANIFEST)

    with patch("httpx.Client") as mock_client_cls:
        client_instance = MagicMock()
        mock_client_cls.return_value.__enter__.return_value = client_instance
        client_instance.get.side_effect = [release_resp, manifest_resp]
        checker.run()

    assert emitted == []


def test_update_checker_silent_on_404():
    """No crash or signal when repo has no releases (private or empty)."""
    from gui.updater import UpdateChecker

    failed = []
    checker = UpdateChecker()
    checker.update_available = MagicMock()
    checker.update_available.emit = lambda _: None
    checker.check_failed = MagicMock()
    checker.check_failed.emit = lambda e: failed.append(e)

    resp_404 = _mock_http_response({}, status=404)
    resp_404.raise_for_status = MagicMock()  # 404 is caught before raise_for_status

    with patch("httpx.Client") as mock_client_cls:
        client_instance = MagicMock()
        mock_client_cls.return_value.__enter__.return_value = client_instance
        client_instance.get.return_value = resp_404
        checker.run()

    assert failed == []


def test_update_checker_emits_check_failed_on_network_error():
    """check_failed signal fires on connection error."""
    from gui.updater import UpdateChecker

    failed = []
    checker = UpdateChecker()
    checker.update_available = MagicMock()
    checker.check_failed = MagicMock()
    checker.check_failed.emit = lambda e: failed.append(e)

    with patch("httpx.Client") as mock_client_cls:
        mock_client_cls.return_value.__enter__.side_effect = ConnectionError("timeout")
        checker.run()

    assert len(failed) == 1
    assert "timeout" in failed[0]


def test_update_checker_no_assets_falls_back_to_html_url():
    """When release has no binary assets, download_url should be the html_url (browser fallback)."""
    from gui.updater import UpdateChecker

    emitted = []
    checker = UpdateChecker()
    checker.update_available = MagicMock()
    checker.update_available.emit = lambda info: emitted.append(info)
    checker.check_failed = MagicMock()

    release_resp = _mock_http_response(MOCK_RELEASE_NO_ASSETS)

    with patch("httpx.Client") as mock_client_cls:
        client_instance = MagicMock()
        mock_client_cls.return_value.__enter__.return_value = client_instance
        client_instance.get.return_value = release_resp
        checker.run()

    assert len(emitted) == 1
    info = emitted[0]
    assert info.can_patch is False
    # download_url must be the html release page, not an empty string
    assert info.download_url == MOCK_RELEASE_NO_ASSETS["html_url"]


# ── DownloadWorker progress ───────────────────────────────────────────────────

def test_download_worker_progress_with_content_length():
    """Progress signal fires correctly when content-length header is present."""
    from gui.updater import DownloadWorker

    progress_vals = []
    worker = DownloadWorker("http://x/file.tar.gz", "file.tar.gz")
    worker.progress = MagicMock()
    worker.progress.emit = lambda v: progress_vals.append(v)
    worker.finished = MagicMock()
    worker.error = MagicMock()

    fake_data = b"A" * 1000

    mock_resp = MagicMock()
    mock_resp.headers = {"content-length": "1000"}
    mock_resp.raise_for_status = MagicMock()
    mock_resp.iter_bytes.return_value = [fake_data[:500], fake_data[500:]]
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)

    with patch("httpx.Client") as mock_client_cls, \
         patch("builtins.open", MagicMock(return_value=MagicMock(
             __enter__=lambda s: MagicMock(write=MagicMock()),
             __exit__=MagicMock(return_value=False)))):
        client_instance = MagicMock()
        mock_client_cls.return_value.__enter__.return_value = client_instance
        client_instance.stream.return_value = mock_resp

        worker.run()

    assert len(progress_vals) >= 1
    assert progress_vals[-1] == 100


def test_download_worker_emits_indeterminate_when_no_content_length():
    """When content-length is absent, progress emits -1 to signal indeterminate mode."""
    from gui.updater import DownloadWorker

    progress_vals = []
    worker = DownloadWorker("http://x/file.tar.gz", "file.tar.gz")
    worker.progress = MagicMock()
    worker.progress.emit = lambda v: progress_vals.append(v)
    worker.finished = MagicMock()
    worker.error = MagicMock()

    fake_data = b"A" * 1000
    mock_resp = MagicMock()
    mock_resp.headers = {}  # no content-length
    mock_resp.raise_for_status = MagicMock()
    mock_resp.iter_bytes.return_value = [fake_data]
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)

    with patch("httpx.Client") as mock_client_cls, \
         patch("builtins.open", MagicMock(return_value=MagicMock(
             __enter__=lambda s: MagicMock(write=MagicMock()),
             __exit__=MagicMock(return_value=False)))):
        client_instance = MagicMock()
        mock_client_cls.return_value.__enter__.return_value = client_instance
        client_instance.stream.return_value = mock_resp

        worker.run()

    assert -1 in progress_vals, "Should emit -1 to trigger indeterminate progress bar"


# ── apply_code_patch safety ───────────────────────────────────────────────────

def test_apply_code_patch_fails_outside_app_bundle():
    """Raises RuntimeError when not running from a .app bundle (expected in dev)."""
    with patch("gui.updater._app_bundle_path", return_value=None):
        with pytest.raises(RuntimeError, match=".app bundle"):
            apply_code_patch("/tmp/fake.tar.gz")


def test_apply_code_patch_path_traversal_blocked():
    """Tar members with ../ paths must not write outside the .app bundle."""
    tmpdir = Path(tempfile.mkdtemp())
    fake_app = tmpdir / "Oligolia.app"
    fake_app.mkdir()

    # Build a malicious tar with a path traversal member
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        info = tarfile.TarInfo(name="../../evil.txt")
        data = b"pwned"
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    buf.seek(0)

    with tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False) as f:
        f.write(buf.read())
        patch_path = f.name

    # The traversal target resolves to tmpdir.parent/evil.txt — pre-clean any leftovers
    evil = tmpdir.parent / "evil.txt"
    if evil.exists():
        evil.unlink()

    try:
        with patch("gui.updater._app_bundle_path", return_value=fake_app), \
             patch("subprocess.run"):
            with pytest.raises(RuntimeError, match="unsafe path"):
                apply_code_patch(patch_path)

        assert not evil.exists(), "Path traversal wrote file outside .app despite fix"
    finally:
        os.unlink(patch_path)
        if evil.exists():
            evil.unlink()


# ── restart_app ───────────────────────────────────────────────────────────────

def test_restart_app_calls_execv():
    """restart_app should call os.execv with the current executable."""
    with patch("os.execv") as mock_execv:
        import sys
        original_argv = sys.argv[:]
        sys.argv = ["oligolia", "--some-flag"]
        try:
            restart_app()
        except Exception:
            pass  # execv replaces process — mock may raise
        sys.argv = original_argv

    mock_execv.assert_called_once()
    args = mock_execv.call_args[0]
    assert args[0] == sys.executable
    assert sys.executable in args[1]
