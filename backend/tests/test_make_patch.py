"""
Tests for scripts/make_patch.py — the producer side of the update system.

The updater rejects a downloaded patch whose sha256 doesn't match the manifest,
so make_patch must record the *actual* checksum of the tarball it just wrote.
This runs the real script against a minimal fake .app bundle via subprocess.
"""

import hashlib
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent


def test_make_patch_records_matching_sha256(tmp_path):
    # Minimal fake .app carrying the core patch files make_patch looks for
    app = tmp_path / "Oligolia.app"
    (app / "Contents" / "MacOS").mkdir(parents=True)
    (app / "Contents" / "Resources").mkdir(parents=True)
    (app / "Contents" / "MacOS" / "Oligolia").write_bytes(b"fake binary payload")
    (app / "Contents" / "Resources" / "version.py").write_text("VERSION = '0.0.0'\n")
    (app / "Contents" / "Info.plist").write_text("<plist></plist>")

    out = tmp_path / "dist"
    script = REPO_ROOT / "scripts" / "make_patch.py"
    result = subprocess.run(
        [sys.executable, str(script), str(app), str(out)],
        cwd=str(REPO_ROOT), capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr

    manifests = list(out.glob("*-manifest.json"))
    patches = list(out.glob("*-mac-patch.tar.gz"))
    assert len(manifests) == 1, "expected exactly one manifest"
    assert len(patches) == 1, "expected exactly one patch tarball"

    manifest = json.loads(manifests[0].read_text())
    assert "patch_sha256" in manifest, "manifest must carry a patch checksum"

    expected = hashlib.sha256(patches[0].read_bytes()).hexdigest()
    assert manifest["patch_sha256"] == expected, \
        "manifest sha256 must match the bytes of the patch it describes"
    assert len(manifest["patch_sha256"]) == 64
