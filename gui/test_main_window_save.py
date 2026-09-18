"""Saving to a place you cannot write must report, not abort.

#55, confirmed path 1: `_save_fasta` / `_save_genbank` called `open(path, "w")`
with no handler, so a PermissionError went straight out of the slot and PyQt6
ended the process with exit -6. The user lost their session for choosing the
wrong folder in a file dialog.

Builds the real MainWindow (~1.7 s offscreen) rather than a stub, because the
thing being checked is the method the menu actions call.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from gui import main_window as mw  # noqa: E402


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def window(app: QApplication) -> "mw.MainWindow":
    return mw.MainWindow()


@pytest.fixture()
def recorded(monkeypatch) -> list[tuple[str, str]]:
    seen: list[tuple[str, str]] = []

    class FakeBox:
        @staticmethod
        def critical(_parent, title, text, *a, **k):
            seen.append((title, text))

        @staticmethod
        def warning(_parent, title, text, *a, **k):
            seen.append((title, text))

        @staticmethod
        def information(_parent, title, text, *a, **k):
            seen.append((title, text))

    monkeypatch.setattr(mw, "QMessageBox", FakeBox)
    return seen


def test_a_good_path_is_written(window, tmp_path, recorded) -> None:
    out = tmp_path / "ok.fasta"
    assert window._write_text(str(out), ">a\nACGT\n") is True
    assert out.read_text() == ">a\nACGT\n"
    assert recorded == [], recorded


def test_an_unwritable_directory_is_reported_not_fatal(window, tmp_path, recorded) -> None:
    """The exact #55 repro, without needing a read-only mount: a missing parent."""
    target = tmp_path / "no-such-dir" / "out.fasta"
    assert window._write_text(str(target), ">a\nACGT\n") is False
    assert len(recorded) == 1, recorded
    title, text = recorded[0]
    assert "save" in title.lower()
    assert str(target) in text
    # It should say what to do about it, not just name the errno.
    assert "different location" in text


def test_a_read_only_directory_is_reported(window, tmp_path, recorded) -> None:
    ro = tmp_path / "ro"
    ro.mkdir()
    ro.chmod(0o500)
    try:
        assert window._write_text(str(ro / "out.fasta"), "x") is False
        assert len(recorded) == 1, recorded
    finally:
        ro.chmod(0o700)  # so tmp_path cleanup works


def test_a_directory_given_as_the_target_is_reported(window, tmp_path, recorded) -> None:
    assert window._write_text(str(tmp_path), "x") is False
    assert len(recorded) == 1, recorded


def test_the_save_actions_route_through_the_guard() -> None:
    """Pin the wiring: no bare open(path, "w") left in the save actions."""
    import ast
    from pathlib import Path

    source = (Path(mw.__file__)).read_text()
    tree = ast.parse(source)
    for fname in ("_save_fasta", "_save_genbank"):
        fn = next(
            (n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == fname),
            None,
        )
        assert fn is not None, f"{fname} not found"
        opens = [
            n for n in ast.walk(fn)
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "open"
        ]
        assert not opens, f"{fname} still calls open() directly; use _write_text (#55)"
