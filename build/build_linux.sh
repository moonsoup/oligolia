#!/bin/bash
# Build Oligolia AppImage for Linux
# Usage: bash build/build_linux.sh
# Output: dist/Oligolia-0.1.0-x86_64.AppImage

set -e
cd "$(dirname "$0")/.."       # project root

VENV="backend/.venv"
PYTHON="$VENV/bin/python"
PYI="$VENV/bin/pyinstaller"
VERSION=$(python3 -c "exec(open('version.py').read()); print(VERSION)")
APP_NAME="Oligolia"
ARCH=$(uname -m)

echo "=== Oligolia Linux build ($ARCH) ==="

# ── Dependencies ────────────────────────────────────────────────────────────
if [ ! -f "$PYTHON" ]; then
    echo "→ Creating virtual environment…"
    python3 -m venv "$VENV"
    "$VENV/bin/pip" install --quiet -r backend/requirements.txt PyInstaller PyQt6 Pillow
fi

# ── 1. Generate icons ───────────────────────────────────────────────────────
echo "→ Generating icons…"
"$PYTHON" assets/make_icon.py

# ── 2. PyInstaller ──────────────────────────────────────────────────────────
echo "→ Building executable…"
"$PYI" build/oligolia.spec \
    --distpath dist \
    --workpath build/work \
    --noconfirm \
    --clean

EXE="dist/${APP_NAME}/${APP_NAME}"
if [ ! -f "$EXE" ]; then
    echo "ERROR: $EXE not found. Check PyInstaller output above."
    exit 1
fi
echo "→ Built: $EXE"

# ── 3. Create .desktop file ─────────────────────────────────────────────────
DESKTOP_FILE="dist/${APP_NAME}/${APP_NAME}.desktop"
cat > "$DESKTOP_FILE" << EOF
[Desktop Entry]
Name=Oligolia
Comment=Gene Editing and Viewing Platform
Exec=Oligolia %f
Icon=oligolia
Terminal=false
Type=Application
Categories=Science;Biology;Education;
MimeType=application/x-fasta;application/x-genbank;text/x-vcf;
Keywords=gene;DNA;RNA;bioinformatics;CRISPR;sequence;
StartupWMClass=Oligolia
EOF

# Copy icon
cp assets/icon.png "dist/${APP_NAME}/oligolia.png"

# ── 4. Create AppImage ───────────────────────────────────────────────────────
if command -v appimagetool &>/dev/null; then
    echo "→ Creating AppImage…"
    APPDIR="build/work/${APP_NAME}.AppDir"
    mkdir -p "$APPDIR/usr/bin" "$APPDIR/usr/share/applications" "$APPDIR/usr/share/icons/hicolor/256x256/apps"

    cp -r "dist/${APP_NAME}/"* "$APPDIR/usr/bin/"
    cp "$DESKTOP_FILE" "$APPDIR/usr/share/applications/"
    cp "assets/icon.png" "$APPDIR/usr/share/icons/hicolor/256x256/apps/oligolia.png"
    cp "assets/icon.png" "$APPDIR/oligolia.png"
    cp "$DESKTOP_FILE" "$APPDIR/"

    APPIMAGE_PATH="dist/${APP_NAME}-${VERSION}-${ARCH}.AppImage"
    ARCH="$ARCH" appimagetool "$APPDIR" "$APPIMAGE_PATH"
    chmod +x "$APPIMAGE_PATH"
    echo ""
    echo "✅  AppImage: $APPIMAGE_PATH"
    echo ""
    echo "User instructions:"
    echo "  1. chmod +x ${APP_NAME}-${VERSION}-${ARCH}.AppImage"
    echo "  2. Double-click or run ./${APP_NAME}-${VERSION}-${ARCH}.AppImage"
    echo "  Optional: right-click → 'Integrate' to add to app menu"
else
    echo "→ appimagetool not found — creating .tar.gz instead"
    echo "  (Download appimagetool from https://github.com/AppImage/AppImageKit/releases)"
    TAR_PATH="dist/${APP_NAME}-${VERSION}-linux-${ARCH}.tar.gz"
    tar czf "$TAR_PATH" -C dist "${APP_NAME}"
    echo ""
    echo "✅  Archive: $TAR_PATH"
    echo ""
    echo "User instructions:"
    echo "  1. tar xzf ${APP_NAME}-${VERSION}-linux-${ARCH}.tar.gz"
    echo "  2. cd ${APP_NAME} && ./${APP_NAME}"
    echo ""
    echo "To also install a .desktop shortcut:"
    echo "  cp ${APP_NAME}.desktop ~/.local/share/applications/"
    echo "  cp oligolia.png ~/.local/share/icons/"
fi

# ── 5. Create .deb package (optional) ───────────────────────────────────────
if command -v dpkg-deb &>/dev/null; then
    echo "→ Creating .deb package…"
    DEB_DIR="build/work/${APP_NAME}-deb"
    mkdir -p "$DEB_DIR/DEBIAN"
    mkdir -p "$DEB_DIR/usr/bin"
    mkdir -p "$DEB_DIR/usr/share/applications"
    mkdir -p "$DEB_DIR/usr/share/icons/hicolor/256x256/apps"
    mkdir -p "$DEB_DIR/usr/lib/${APP_NAME,,}"

    # Control file
    cat > "$DEB_DIR/DEBIAN/control" << EOF
Package: oligolia
Version: ${VERSION}
Architecture: ${ARCH}
Maintainer: Oligolia Project <support@oligolia.app>
Description: Gene Editing and Viewing Platform
 Advanced bioinformatics tool for sequence editing, CRISPR design,
 database search (NCBI, Ensembl, UniProt, KEGG), alignment, and more.
Depends: libgl1
EOF

    cp -r "dist/${APP_NAME}/"* "$DEB_DIR/usr/lib/${APP_NAME,,}/"
    ln -sf "/usr/lib/${APP_NAME,,}/${APP_NAME}" "$DEB_DIR/usr/bin/oligolia"
    cp "$DESKTOP_FILE" "$DEB_DIR/usr/share/applications/"
    cp "assets/icon.png" "$DEB_DIR/usr/share/icons/hicolor/256x256/apps/oligolia.png"

    DEB_PATH="dist/${APP_NAME,,}_${VERSION}_${ARCH}.deb"
    dpkg-deb --build "$DEB_DIR" "$DEB_PATH"
    echo "✅  .deb: $DEB_PATH"
fi
