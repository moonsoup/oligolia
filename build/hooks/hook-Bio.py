"""PyInstaller hook for Biopython — collect all submodules and data."""
from PyInstaller.utils.hooks import collect_all

datas, binaries, hiddenimports = collect_all("Bio")
