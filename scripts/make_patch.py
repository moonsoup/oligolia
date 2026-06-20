"""
Generate a code-only patch from a built .app bundle.

A "code patch" contains only the files that change between releases:
  - Contents/MacOS/Oligolia        (11 MB — compiled Python code archive)
  - Contents/Resources/version.py  (4 KB — version number)
  - Contents/Info.plist            (version strings)

This is ~11 MB vs ~260 MB for a full DMG — 96% smaller.

The Frameworks/ directory (Qt + Python dylibs, ~600 MB) is excluded
because it only changes when PyQt6 or Python itself is upgraded.

Usage:
    python scripts/make_patch.py [app_path] [output_dir]

Defaults:
    app_path   = dist/Oligolia.app
    output_dir = dist/
"""

import sys
import tarfile
import json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from version import VERSION

APP_PATH = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("dist/Oligolia.app")
OUT_DIR = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("dist")

# Files included in a code patch (relative to .app root)
PATCH_FILES = [
    "Contents/MacOS/Oligolia",
    "Contents/Resources/version.py",
    "Contents/Info.plist",
]

# Optional: include any changed Resources (not stdlib, not frameworks)
PATCH_RESOURCES_PATTERNS = [
    "Contents/Resources/backend",
    "Contents/Resources/gui",
    "Contents/Resources/oligolia.py",
]


def make_patch() -> None:
    if not APP_PATH.exists():
        print(f"ERROR: {APP_PATH} not found. Run 'make mac' first.")
        sys.exit(1)

    OUT_DIR.mkdir(exist_ok=True)
    patch_name = f"Oligolia-{VERSION}-mac-patch.tar.gz"
    patch_path = OUT_DIR / patch_name
    manifest_path = OUT_DIR / f"Oligolia-{VERSION}-manifest.json"

    print(f"Creating code patch: {patch_path}")

    included: list[Path] = []
    total_bytes = 0

    with tarfile.open(patch_path, "w:gz", compresslevel=9) as tar:
        # Core patch files
        for rel in PATCH_FILES:
            p = APP_PATH / rel
            if p.exists():
                tar.add(p, arcname=rel)
                sz = p.stat().st_size
                total_bytes += sz
                included.append(p)
                print(f"  + {rel}  ({sz / 1024 / 1024:.1f} MB)")

        # Our Python source directories (if present as .py files)
        for pattern in PATCH_RESOURCES_PATTERNS:
            p = APP_PATH / pattern
            if p.exists():
                tar.add(p, arcname=pattern.removeprefix("Contents/"))
                sz = sum(f.stat().st_size for f in p.rglob("*") if f.is_file()) if p.is_dir() else p.stat().st_size
                total_bytes += sz
                print(f"  + {pattern}  ({sz / 1024:.0f} KB)")

    patch_size = patch_path.stat().st_size
    print(f"\nPatch size:  {patch_size / 1024 / 1024:.1f} MB")

    # Generate manifest
    manifest = {
        "version": VERSION,
        "requires_full": False,
        "min_compatible_base": "0.3.0",  # oldest version this patch applies to
        "changelog": f"Oligolia {VERSION}",
        "assets": {
            "darwin_patch": patch_name,
            "darwin_full": f"Oligolia-{VERSION}-mac.dmg",
            "windows_full": f"Oligolia-{VERSION}-Setup.exe",
            "linux_full": f"Oligolia-{VERSION}-x86_64.AppImage",
        },
    }
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"Manifest:    {manifest_path}")
    print("\nSavings vs full DMG: ~96% smaller")


if __name__ == "__main__":
    make_patch()
