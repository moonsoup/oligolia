"""PyInstaller hook for PyQt6 — ensure all needed Qt plugins are bundled."""
from PyInstaller.utils.hooks import collect_all

datas, binaries, hiddenimports = collect_all("PyQt6")
