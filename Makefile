# Oligolia build targets
# Usage: make mac | make win | make linux | make all

PYTHON  := backend/.venv/bin/python
PYI     := backend/.venv/bin/pyinstaller
VERSION := $(shell $(PYTHON) -c "from version import VERSION; print(VERSION)")

SPEC    := build/oligolia.spec
DIST    := dist
WORK    := build/work

.PHONY: all mac win linux clean icons venv test

# ── Default: local platform ────────────────────────────────────────────────
all: icons
	@case "$$(uname -s)" in \
	  Darwin*)  $(MAKE) mac ;; \
	  Linux*)   $(MAKE) linux ;; \
	  *)        $(MAKE) win ;; \
	esac

# ── macOS DMG + code patch ─────────────────────────────────────────────────
mac: icons _pyinstaller
	@echo "→ Creating DMG…"
	@TMP=$$(mktemp -d) && \
	cp -R $(DIST)/Oligolia.app $$TMP/ && \
	ln -s /Applications $$TMP/Applications && \
	hdiutil create \
	    -volname "Oligolia $(VERSION)" \
	    -srcfolder $$TMP \
	    -ov -format UDZO -imagekey zlib-level=9 \
	    $(DIST)/Oligolia-$(VERSION)-mac.dmg 2>&1 | tail -1 && \
	rm -rf $$TMP
	@echo "✅  $(DIST)/Oligolia-$(VERSION)-mac.dmg"
	@echo "→ Creating code patch…"
	$(PYTHON) scripts/make_patch.py $(DIST)/Oligolia.app $(DIST)
	@echo "✅  $(DIST)/Oligolia-$(VERSION)-mac-patch.tar.gz"

# ── Windows (run on Windows or in CI) ─────────────────────────────────────
win: icons _pyinstaller
	@echo "→ Creating Windows installer (requires Inno Setup)…"
	@iscc build/inno_setup.iss 2>/dev/null || \
	  (powershell -Command "Compress-Archive -Path '$(DIST)/Oligolia/*' \
	    -DestinationPath '$(DIST)/Oligolia-$(VERSION)-win.zip' -Force" && \
	   echo "✅  $(DIST)/Oligolia-$(VERSION)-win.zip (zip fallback)")

# ── Linux AppImage ─────────────────────────────────────────────────────────
linux: icons _pyinstaller
	@echo "→ Creating AppImage…"
	@bash build/build_linux.sh

# ── Shared PyInstaller step ────────────────────────────────────────────────
_pyinstaller:
	@echo "→ Running PyInstaller (version $(VERSION))…"
	$(PYI) $(SPEC) --distpath $(DIST) --workpath $(WORK) --noconfirm --clean

# ── Icons ──────────────────────────────────────────────────────────────────
icons:
	@echo "→ Generating icons…"
	$(PYTHON) assets/make_icon.py

# ── Virtual environment ────────────────────────────────────────────────────
venv:
	python3 -m venv backend/.venv
	backend/.venv/bin/pip install -r backend/requirements.txt PyInstaller PyQt6 Pillow

# ── Tests ──────────────────────────────────────────────────────────────────
test:
	cd backend && .venv/bin/python -m pytest tests/ -q

# ── Clean build artifacts ──────────────────────────────────────────────────
clean:
	rm -rf $(DIST) $(WORK)
	@echo "Cleaned dist/ and build/work/"
