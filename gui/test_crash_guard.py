"""#55: an unhandled exception in a slot must not kill the process.

PyQt6 aborts (exit -6) when an exception escapes a slot, so before this every
missed edge case was a crash rather than an error dialog. Two confirmed paths:
saving to an unwritable location, and `[]` as workflow step parameters.

These tests drive the hook directly rather than through Qt, because the thing
under test is what happens when Qt is already unwinding.
"""

from __future__ import annotations

import sys

import pytest

from gui import crash_guard


@pytest.fixture(autouse=True)
def _forget_seen_faults():
    """The dedup guard keeps module state, so reset it around every test.

    Without this the suite is order-dependent: whichever test raises a given
    fault first gets the dialog and the rest silently see none.
    """
    crash_guard.reset()
    yield
    crash_guard.reset()


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


# --- the live Qt question source inspection cannot answer ---
#
# Codex, 2026-09-17: "source inspection cannot establish whether the installed
# PyQt6 version routes an exception escaping a slot through sys.excepthook before
# aborting", and "no re-entrancy guard or exception deduplication exists".
#
# Both are settled here by running a real QApplication in a subprocess. The first
# question needed an answer because the #55 commit claimed the abort is prevented;
# the second turned out to be a real defect — a QTimer repeating every 5 ms drove
# the hook 55 times in 400 ms, i.e. 55 modal dialogs.

_PROBE = '''
import os, sys
os.environ["QT_QPA_PLATFORM"] = "offscreen"
sys.path.insert(0, {root!r})

from PyQt6.QtCore import QTimer
from PyQt6.QtWidgets import QApplication
from gui.crash_guard import install

app = QApplication([])
seen = []
install(show=lambda title, body: seen.append(body))

def boom():
    raise PermissionError(13, "Permission denied", "/nope/out.fasta")

def after():
    print("RESULT", len(seen))
    app.quit()

{schedule}
QTimer.singleShot(400, after)
sys.exit(app.exec())
'''


def _run_probe(tmp_path, schedule: str) -> tuple[int, int]:
    """Run a real Qt app that raises from a slot. Returns (exit code, dialogs)."""
    import subprocess
    import sys as _sys
    from pathlib import Path

    root = str(Path(__file__).resolve().parents[1])
    script = tmp_path / "probe.py"
    script.write_text(_PROBE.format(root=root, schedule=schedule))

    done = subprocess.run(
        [_sys.executable, str(script)], capture_output=True, text=True, timeout=120
    )
    count = -1
    for line in done.stdout.splitlines():
        if line.startswith("RESULT"):
            count = int(line.split()[1])
    return done.returncode, count


def test_an_exception_in_a_real_slot_reaches_the_hook_and_the_app_lives(tmp_path) -> None:
    """#55's central claim, executed rather than asserted.

    Without the hook PyQt6 ends the process with SIGABRT (-6). If this test ever
    fails with a negative return code, the hook is not intercepting and #55 should
    be reopened.
    """
    code, dialogs = _run_probe(tmp_path, "QTimer.singleShot(10, boom)")
    assert code == 0, f"process did not survive the slot exception (exit {code})"
    assert dialogs == 1, dialogs


def test_a_slot_that_raises_every_tick_shows_one_dialog_not_fifty(tmp_path) -> None:
    """The defect Codex predicted: measured at 55 dialogs in 400 ms before the guard."""
    code, dialogs = _run_probe(
        tmp_path, "t = QTimer(); t.timeout.connect(boom); t.start(5)"
    )
    assert code == 0, f"process did not survive the repeating slot (exit {code})"
    assert dialogs == 1, f"{dialogs} dialogs for one repeating fault — that is a loop"


# --- the dedup guard, driven directly ---

def test_the_same_fault_twice_shows_one_dialog(tmp_path) -> None:
    crash_guard.reset()
    shown: list[str] = []
    e = _boom()
    for _ in range(10):
        crash_guard.handle_exception(
            type(e), e, e.__traceback__,
            show=lambda t, b: shown.append(b), log_path=tmp_path / "c.log",
        )
    assert len(shown) == 1, len(shown)


def test_a_different_fault_still_gets_its_own_dialog(tmp_path) -> None:
    crash_guard.reset()
    shown: list[str] = []
    log = tmp_path / "c.log"

    e1 = _boom()
    crash_guard.handle_exception(type(e1), e1, e1.__traceback__,
                                 show=lambda t, b: shown.append(b), log_path=log)
    try:
        raise ValueError("something else entirely")
    except ValueError as e2:
        crash_guard.handle_exception(type(e2), e2, e2.__traceback__,
                                     show=lambda t, b: shown.append(b), log_path=log)

    assert len(shown) == 2, shown
    assert "PermissionError" in shown[0]
    assert "ValueError" in shown[1]


def test_a_repeated_fault_stops_filling_the_log(tmp_path) -> None:
    crash_guard.reset()
    log = tmp_path / "c.log"
    e = _boom()
    for _ in range(50):
        crash_guard.handle_exception(type(e), e, e.__traceback__,
                                     show=lambda t, b: None, log_path=log)
    assert log.read_text().count("Traceback") == crash_guard.MAX_LOGS_PER_FAULT
    assert "will not be logged" in log.read_text()
