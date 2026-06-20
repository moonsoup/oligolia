#!/bin/bash
# release.sh — bump version, build locally, tag, push → GitHub Actions builds all platforms
#
# Usage:
#   bash scripts/release.sh 0.2.0          # release version 0.2.0
#   bash scripts/release.sh 0.2.0 --dry-run  # show what would happen, don't push

set -e
cd "$(dirname "$0")/.."   # project root

NEW_VERSION="${1:-}"
DRY_RUN="${2:-}"

if [ -z "$NEW_VERSION" ]; then
    echo "Usage: bash scripts/release.sh <version>  [--dry-run]"
    echo "Example: bash scripts/release.sh 0.2.0"
    exit 1
fi

PYTHON="backend/.venv/bin/python"
CURRENT=$("$PYTHON" -c "from version import VERSION; print(VERSION)" 2>/dev/null || echo "unknown")

echo "=== Oligolia Release Script ==="
echo "Current version : $CURRENT"
echo "New version     : $NEW_VERSION"
echo ""

if [ "$DRY_RUN" = "--dry-run" ]; then
    echo "[DRY RUN] Would:"
    echo "  1. Update version.py  → $NEW_VERSION"
    echo "  2. make mac            → dist/Oligolia-$NEW_VERSION-mac.dmg"
    echo "  3. git commit + tag   → v$NEW_VERSION"
    echo "  4. git push + push tag → triggers GitHub Actions (Win + Linux builds)"
    exit 0
fi

# ── 1. Bump version.py (single source of truth — spec/ISS/package.json all read from here) ──
echo "→ Bumping version.py  $CURRENT → $NEW_VERSION…"
sed -i '' "s/VERSION = \".*\"/VERSION = \"$NEW_VERSION\"/" version.py

# Keep frontend package.json in sync (cosmetic, doesn't affect the build)
sed -i '' "s/\"version\": \".*\"/\"version\": \"$NEW_VERSION\"/" frontend/package.json

# Keep Inno Setup in sync (it doesn't read version.py)
sed -i '' "s/#define MyAppVersion \".*\"/#define MyAppVersion \"$NEW_VERSION\"/" build/inno_setup.iss

echo "  version.py     → $(grep 'VERSION = ' version.py)"

# ── 2. Build locally (macOS DMG) ────────────────────────────────────────────
echo "→ Building macOS app + DMG…"
make mac

DMG="dist/Oligolia-$NEW_VERSION-mac.dmg"
if [ ! -f "$DMG" ]; then
    echo "ERROR: DMG not found at $DMG"
    exit 1
fi
echo "→ Built: $DMG ($(du -sh "$DMG" | cut -f1))"

# ── 3. Commit + tag ─────────────────────────────────────────────────────────
echo "→ Committing version bump…"
git add version.py frontend/package.json build/inno_setup.iss
git commit -m "Release v$NEW_VERSION

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"

echo "→ Tagging v$NEW_VERSION…"
git tag "v$NEW_VERSION"

# ── 4. Push ─────────────────────────────────────────────────────────────────
echo "→ Pushing to GitHub…"
git push origin main
git push origin "v$NEW_VERSION"

echo ""
echo "✅  Release v$NEW_VERSION pushed!"
echo ""
echo "GitHub Actions is now building:"
echo "  • Windows installer  (Oligolia-$NEW_VERSION-Setup.exe)"
echo "  • Linux AppImage     (Oligolia-$NEW_VERSION-x86_64.AppImage)"
echo ""
echo "macOS DMG is already built locally at:"
echo "  $DMG"
echo ""
echo "Watch the build: https://github.com/moonsoup/oligolia/actions"
echo "Release page:    https://github.com/moonsoup/oligolia/releases/tag/v$NEW_VERSION"
