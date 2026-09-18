"""#55: an unhandled exception in a slot must not kill the process.

PyQt6 aborts (exit -6) when an exception escapes a slot, so before this every
missed edge case was a crash rather than an error dialog. Two confirmed paths:
saving to an unwritable location, and `[]` as workflow step parameters.

These tests drive the hook directly rather than through Qt, because the thing
under test is what happens when Qt is already unwinding.
"""

from __future__ import annotations

import sys

from gui import crash_guard


def _boom() -> Exception:
    try:
        raise PermissionError(13, "Permission denied", "/nope/out.fasta")
    except PermissionError as e:
        return e


def test_the_hook_shows_the_error_instead_of_dying(tmp_path) -> None:
    shown: list[tuple[str, str]] = []
    log = tmp_path / "crash.log"
    e = _boom()

    crash_guard.handle_exception(
        type(e), e, e.__traceback__,
        show=lambda title, body: shown.append((title, body)),
        log_path=log,
    )

    assert len(shown) == 1, shown
    title, body = shown[0]
    assert "PermissionError" in body
    assert "Permission denied" in body
    # The user needs to know where the detail went.
    assert str(log) in body


def test_the_full_traceback_goes_to_the_log(tmp_path) -> None:
    log = tmp_path / "crash.log"
    e = _boom()

    crash_guard.handle_exception(type(e), e, e.__traceback__, show=lambda t, b: None, log_path=log)

    text = log.read_text()
    assert "Traceback" in text
    assert "PermissionError" in text
    assert "_boom" in text, "the frame that raised should be in the traceback"


def test_the_log_appends_rather_than_truncating(tmp_path) -> None:
    log = tmp_path / "crash.log"
    e = _boom()
    for _ in range(3):
        crash_guard.handle_exception(type(e), e, e.__traceback__, show=lambda t, b: None, log_path=log)
    assert log.read_text().count("Traceback") == 3


def test_keyboard_interrupt_is_left_to_the_default_hook(tmp_path) -> None:
    """Ctrl-C must still stop the program; swallowing it would be worse than crashing."""
    shown: list[tuple[str, str]] = []
    delegated: list[bool] = []
    e = KeyboardInterrupt()

    crash_guard.handle_exception(
        KeyboardInterrupt, e, None,
        show=lambda t, b: shown.append((t, b)),
        log_path=tmp_path / "crash.log",
        fallback=lambda *a: delegated.append(True),
    )

    assert shown == [], "Ctrl-C should not raise a dialog"
    assert delegated == [True]


def test_a_failure_inside_the_dialog_does_not_propagate(tmp_path) -> None:
    """The hook runs while the interpreter is already unwinding; it must not raise."""
    e = _boom()
    crash_guard.handle_exception(
        type(e), e, e.__traceback__,
        show=lambda t, b: (_ for _ in ()).throw(RuntimeError("no display")),
        log_path=tmp_path / "crash.log",
    )  # must return normally


def test_an_unwritable_log_does_not_propagate(tmp_path) -> None:
    """Logging to a bad path must not turn the crash guard into the crash."""
    e = _boom()
    crash_guard.handle_exception(
        type(e), e, e.__traceback__,
        show=lambda t, b: None,
        log_path=tmp_path / "no-such-dir" / "deeper" / "crash.log",
    )


def test_install_sets_and_restores_the_hook() -> None:
    before = sys.excepthook
    restore = crash_guard.install(show=lambda t, b: None)
    try:
        assert sys.excepthook is not before
    finally:
        restore()
    assert sys.excepthook is before


def test_the_entry_point_installs_the_hook() -> None:
    """The guard is worthless if nothing turns it on.

    Resolved through the AST so the import may be aliased however reads best.
    """
    import ast
    from pathlib import Path

    source = (Path(__file__).resolve().parents[1] / "oligolia.py").read_text()
    tree = ast.parse(source)

    # Names bound by `from gui.crash_guard import install [as X]`.
    bound = {
        (alias.asname or alias.name)
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and (node.module or "").endswith("crash_guard")
        for alias in node.names
        if alias.name == "install"
    }
    assert bound, "oligolia.py must import install from gui.crash_guard (#55)"

    called = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    assert bound & called, f"imported {bound} but never called it (#55)"
