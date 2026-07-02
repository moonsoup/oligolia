# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller spec for Oligolia — Gene Editing Platform.

Usage (run from project root):
    pyinstaller build/oligolia.spec
"""

import sys
import os
from pathlib import Path

project_root = str(Path(SPECPATH).parent)
sys.path.insert(0, project_root)
from version import VERSION  # single source of truth

block_cipher = None

# ── Hidden imports ────────────────────────────────────────────────────────────
# Modules PyInstaller can't detect automatically (dynamic imports, plugins)

HIDDEN_IMPORTS = [
    # PyQt6
    "PyQt6",
    "PyQt6.QtCore",
    "PyQt6.QtGui",
    "PyQt6.QtWidgets",
    "PyQt6.sip",
    # Biopython — many modules are imported dynamically
    "Bio",
    "Bio.Seq",
    "Bio.SeqIO",
    "Bio.SeqRecord",
    "Bio.SeqFeature",
    "Bio.Align",
    "Bio.Align.substitution_matrices",
    "Bio.GenBank",
    "Bio.GenBank.Scanner",
    "Bio.File",
    "Bio.Data",
    "Bio.Data.CodonTable",
    "Bio.Data.IUPACData",
    "Bio.Alphabet",
    # Version comparison
    "packaging",
    "packaging.version",
    # Network
    "httpx",
    "httpx._transports",
    "httpx._transports.default",
    "httpcore",
    "anyio",
    "sniffio",
    # Pydantic
    "pydantic",
    "pydantic.v1",
    "pydantic_core",
    # Stdlib extras
    "email",
    "email.mime",
    "email.mime.multipart",
    "email.mime.text",
    "email.mime.base",
    # Backend and GUI packages
    "backend",
    "backend.models",
    "backend.services",
    "backend.formats",
    "backend.routers",
    "gui",
    "gui.panels",
    "gui.workers",
    "gui.styles",
]

# ── Data files ────────────────────────────────────────────────────────────────

import Bio as _bio
bio_path = str(Path(_bio.__file__).parent)

datas = [
    # Biopython data (codon tables, substitution matrices, etc.)
    (os.path.join(bio_path, "Data"), "Bio/Data"),
    (os.path.join(bio_path, "Align"), "Bio/Align"),
    # Application assets
    (os.path.join(project_root, "assets", "icon.png"), "assets"),
    # Bundled fonts (issue #46) — loaded at startup via QFontDatabase
    (os.path.join(project_root, "assets", "fonts"), "assets/fonts"),
    # Version file — bundled so the app knows its own version at runtime
    (os.path.join(project_root, "version.py"), "."),
    # Bundled reference library of common vector parts (issue #42)
    (os.path.join(project_root, "backend", "data"), "backend/data"),
]

# ── Analysis ──────────────────────────────────────────────────────────────────

a = Analysis(
    [os.path.join(project_root, "oligolia.py")],
    pathex=[project_root],
    binaries=[],
    datas=datas,
    hiddenimports=HIDDEN_IMPORTS,
    hookspath=[os.path.join(project_root, "build", "hooks")],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        # Exclude heavy unused modules
        "tkinter",
        "matplotlib",
        "scipy",
        "sklearn",
        "tensorflow",
        "torch",
        "notebook",
        "IPython",
        "jupyter",
        "test",
        "unittest",
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# ── Executable ────────────────────────────────────────────────────────────────

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="Oligolia",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,            # No terminal window
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=os.path.join(project_root, "assets", "icon.ico"),
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="Oligolia",
)

# ── macOS .app bundle ─────────────────────────────────────────────────────────

if sys.platform == "darwin":
    app = BUNDLE(
        coll,
        name="Oligolia.app",
        icon=os.path.join(project_root, "assets", "icon.icns"),
        bundle_identifier="com.oligolia.geneplatform",
        version=VERSION,
        info_plist={
            "NSPrincipalClass": "NSApplication",
            "NSAppleScriptEnabled": False,
            "CFBundleDisplayName": "Oligolia",
            "CFBundleName": "Oligolia",
            "CFBundleVersion": "0.1.0",
            "CFBundleShortVersionString": "0.1.0",
            "NSHighResolutionCapable": True,
            "LSApplicationCategoryType": "public.app-category.medical",
            "CFBundleDocumentTypes": [
                {
                    "CFBundleTypeName": "FASTA sequence",
                    "CFBundleTypeExtensions": ["fasta", "fa", "fna", "faa"],
                    "CFBundleTypeRole": "Editor",
                },
                {
                    "CFBundleTypeName": "GenBank sequence",
                    "CFBundleTypeExtensions": ["gb", "gbk"],
                    "CFBundleTypeRole": "Editor",
                },
                {
                    "CFBundleTypeName": "VCF variant file",
                    "CFBundleTypeExtensions": ["vcf"],
                    "CFBundleTypeRole": "Viewer",
                },
            ],
        },
    )
