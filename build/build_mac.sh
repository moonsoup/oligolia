#!/bin/bash
# Build Oligolia.dmg for macOS
# Usage: bash build/build_mac.sh
# Output: dist/Oligolia-0.1.0-mac.dmg

set -e
cd "$(dirname "$0")/.."       # project root

VENV="backend/.venv"
PYTHON="$VENV/bin/python"
PYI="$VENV/bin/pyinstaller"
VERSION="0.1.0"
APP_NAME="Oligolia"
DMG_NAME="${APP_NAME}-${VERSION}-mac.dmg"

echo "=== Oligolia macOS build ==="
echo "Python: $("$PYTHON" --version)"

# ── 1. Generate icons ───────────────────────────────────────────────────────
echo "→ Generating icons…"
"$PYTHON" assets/make_icon.py

# ── 2. PyInstaller ──────────────────────────────────────────────────────────
echo "→ Building .app bundle (this takes 2–5 minutes)…"
"$PYI" build/oligolia.spec \
    --distpath dist \
    --workpath build/work \
    --noconfirm \
    --clean

APP_PATH="dist/${APP_NAME}.app"
if [ ! -d "$APP_PATH" ]; then
    echo "ERROR: $APP_PATH not created. Check PyInstaller output above."
    exit 1
fi
echo "→ Built: $APP_PATH"

# ── 3. Code sign (skip if no Developer ID) ──────────────────────────────────
if command -v codesign &>/dev/null && security find-identity -v -p codesigning 2>/dev/null | grep -q "Developer ID"; then
    IDENTITY=$(security find-identity -v -p codesigning | grep "Developer ID" | head -1 | awk -F'"' '{print $2}')
    echo "→ Code signing with: $IDENTITY"
    codesign --force --deep --sign "$IDENTITY" "$APP_PATH"
else
    echo "→ Skipping code signing (no Developer ID certificate found)"
    echo "  Users may need to right-click → Open on first launch."
fi

# ── 4. Create DMG ─────────────────────────────────────────────────────────
echo "→ Creating DMG…"
DMG_PATH="dist/$DMG_NAME"

# Use hdiutil (no extra tools required)
TMP_DIR="$(mktemp -d)"
cp -R "$APP_PATH" "$TMP_DIR/"
# Create a symlink to /Applications for drag-and-drop install
ln -s /Applications "$TMP_DIR/Applications"

hdiutil create \
    -volname "${APP_NAME} ${VERSION}" \
    -srcfolder "$TMP_DIR" \
    -ov \
    -format UDZO \
    -imagekey zlib-level=9 \
    "$DMG_PATH"

rm -rf "$TMP_DIR"
echo ""
echo "✅  Done! Installer: $DMG_PATH"
echo ""
echo "User instructions:"
echo "  1. Double-click Oligolia-${VERSION}-mac.dmg"
echo "  2. Drag Oligolia.app to Applications"
echo "  3. Launch from Launchpad or Applications folder"
