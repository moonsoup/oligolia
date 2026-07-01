#!/usr/bin/env bash
# Build a 0.3.1 code patch from the existing Oligolia.app bundle.
# The patch updates the main binary and version.py inside the installed .app.
set -e

VER="0.3.1"
DIST="$(cd "$(dirname "$0")/.." && pwd)/dist"
APP="$DIST/Oligolia.app"
PATCH_NAME="Oligolia-${VER}-mac-patch.tar.gz"
MANIFEST_NAME="Oligolia-${VER}-manifest.json"

if [ ! -d "$APP" ]; then
  echo "ERROR: $APP not found — build the app first"
  exit 1
fi

echo "Building patch $PATCH_NAME from $APP ..."

# Update version.py inside the bundle to 0.3.1
VERSION_PY="$APP/Contents/Resources/version.py"
if [ -f "$VERSION_PY" ]; then
  sed -i '' 's/VERSION = "0\.3\.0"/VERSION = "0.3.1"/' "$VERSION_PY"
  echo "  Updated $VERSION_PY"
fi

# Build tar with only the files the patch mechanism expects
cd "$APP/.."
tar -czf "$DIST/$PATCH_NAME" \
  --exclude='*.pyc' \
  Oligolia.app/Contents/MacOS/Oligolia \
  Oligolia.app/Contents/Resources/version.py

# Strip leading "Oligolia.app/" so member.name matches apply_code_patch expectations
# (apply_code_patch does: dest = app / member.name)
# Rebuild with correct paths
cd "$APP"
tar -czf "$DIST/$PATCH_NAME" \
  --exclude='*.pyc' \
  Contents/MacOS/Oligolia \
  Contents/Resources/version.py

echo "  Created $DIST/$PATCH_NAME"

# Write manifest
cat > "$DIST/$MANIFEST_NAME" <<EOF
{
  "version": "${VER}",
  "requires_full": false,
  "min_compatible_base": "0.3.0",
  "changelog": "Oligolia ${VER} — bug fix release",
  "assets": {
    "darwin_patch": "${PATCH_NAME}",
    "darwin_full": "Oligolia-${VER}-mac.dmg",
    "windows_full": "Oligolia-${VER}-Setup.exe",
    "linux_full": "Oligolia-${VER}-x86_64.AppImage"
  }
}
EOF
echo "  Created $DIST/$MANIFEST_NAME"
echo "Done."
