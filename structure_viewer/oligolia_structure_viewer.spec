# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller spec for the Oligolia Structure Viewer companion app.

Deliberately separate from build/oligolia.spec — this is the ONLY place
PyQt6-WebEngine (and its bundled Chromium runtime) gets pulled in, so the core
Oligolia installer never grows a Chromium dependency.

Usage (run from project root):
    pyinstaller structure_viewer/oligolia_structure_viewer.spec
"""

import sys
import os
from pathlib import Path

project_root = str(Path(SPECPATH).parent)
viewer_root = os.path.join(project_root, "structure_viewer")
sys.path.insert(0, viewer_root)
from version import VERSION  # noqa: E402 — structure_viewer/version.py, independent of the core app's

block_cipher = None

HIDDEN_IMPORTS = [
    "PyQt6",
    "PyQt6.QtCore",
    "PyQt6.QtGui",
    "PyQt6.QtWidgets",
    "PyQt6.QtWebEngineWidgets",
    "PyQt6.QtWebEngineCore",
    "PyQt6.QtWebChannel",
    "PyQt6.sip",
]

datas = [
    (os.path.join(viewer_root, "assets", "3Dmol-min.js"), "assets"),
    (os.path.join(viewer_root, "version.py"), "."),
]

a = Analysis(
    [os.path.join(viewer_root, "main.py")],
    pathex=[viewer_root],
    binaries=[],
    datas=datas,
    hiddenimports=HIDDEN_IMPORTS,
    hookspath=[os.path.join(project_root, "build", "hooks")],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "tkinter", "matplotlib", "scipy", "sklearn", "tensorflow", "torch",
        "notebook", "IPython", "jupyter", "test", "unittest",
        "Bio",  # this app has no bioinformatics dependency — just renders a .pdb file
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="Oligolia Structure Viewer",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
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
    name="Oligolia Structure Viewer",
)

if sys.platform == "darwin":
    app = BUNDLE(
        coll,
        name="Oligolia Structure Viewer.app",
        icon=os.path.join(project_root, "assets", "icon.icns"),
        bundle_identifier="com.oligolia.structureviewer",
        version=VERSION,
        info_plist={
            "NSPrincipalClass": "NSApplication",
            "CFBundleDisplayName": "Oligolia Structure Viewer",
            "CFBundleName": "Oligolia Structure Viewer",
            "CFBundleVersion": VERSION,
            "CFBundleShortVersionString": VERSION,
            "NSHighResolutionCapable": True,
            "LSApplicationCategoryType": "public.app-category.medical",
            "LSUIElement": False,
        },
    )
