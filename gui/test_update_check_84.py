"""#84 — building a MainWindow must not start a live update check in tests.

`MainWindow.__init__` unconditionally called `_start_update_check()`, so every GUI
test that built the real window fired an `UpdateChecker` QThread at the GitHub
releases API. Two consequences the Codex verifier on #82 hit: the suite makes an
outbound request (its no-network sandbox logged a handled DNS failure), and — when
the code under test is older than the latest published release — `update_available`
fires and `_on_update_available` calls `dlg.exec()`, a *modal* dialog that nobody
can dismiss under `QT_QPA_PLATFORM=offscreen`.

The acceptance list on the issue is the spec:

1. constructing `MainWindow` in the test suite makes no `UpdateChecker`
2. default construction still starts the check — the opt-out is additive
3. forcing `update_available` with `check_updates=False` shows no dialog
4. (kept alongside) Help → Check for Updates is unaffected by the opt-out

`UpdateChecker.start` is patched in every test here, so no test in this file can
reach the network even if the guard under test regresses, and `UpdateDialog` is
replaced by a recorder so no real modal dialog can ever open.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from gui import main_window as mw  # noqa: E402
from gui.updater import UpdateChecker, UpdateInfo  # noqa: E402

#: The environment opt-out gui/conftest.py sets for the whole GUI session.
ENV = "OLIGOLIA_NO_UPDATE_CHECK"


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def started(monkeypatch) -> list[UpdateChecker]:
    """Records every `UpdateChecker.start()` instead of running the thread."""
    seen: list[UpdateChecker] = []
    monkeypatch.setattr(UpdateChecker, "start", lambda self: seen.append(self))
    return seen


@pytest.fixture()
def dialogs(monkeypatch) -> list[UpdateInfo]:
    """Records every UpdateDialog that would have been opened, and opens none."""
    shown: list[UpdateInfo] = []

    class FakeDialog:
        def __init__(self, info: UpdateInfo, parent=None) -> None:
            self._info = info

        def setStyleSheet(self, _sheet: str) -> None:
            pass

        def exec(self) -> int:
            shown.append(self._info)
            return 0

    monkeypatch.setattr(mw, "UpdateDialog", FakeDialog)
    return shown


def _info() -> UpdateInfo:
    return UpdateInfo(
        version="99.0.0", body="notes", html_url="https://example.invalid/r",
        patch_url="", full_url="https://example.invalid/d",
        requires_full=True, min_compatible_base="0.0.0",
    )


# ── 1. the suite starts no update checker ────────────────────────────────────

def test_suite_construction_starts_no_update_checker(app, started, monkeypatch):
    """Under the GUI session's env opt-out, a plain MainWindow() checks nothing.

    This is the bullet that protects tests nobody has written yet: gui/conftest.py
    sets the variable for the whole session, so a new test that forgets
    `check_updates=False` still makes no request.
    """
    assert os.environ.get(ENV) == "1", "gui/conftest.py must set the session opt-out"
    win = mw.MainWindow()
    assert started == []
    assert win._update_checker is None


def test_explicit_opt_out_starts_no_update_checker(app, started, monkeypatch):
    """`check_updates=False` wins even with no environment variable set."""
    monkeypatch.delenv(ENV, raising=False)
    win = mw.MainWindow(check_updates=False)
    assert started == []
    assert win._update_checker is None


def test_calling_start_update_check_directly_is_a_no_op_when_opted_out(app, started, monkeypatch):
    """The guard lives in `_start_update_check`, not only at the call site."""
    monkeypatch.delenv(ENV, raising=False)
    win = mw.MainWindow(check_updates=False)
    win._start_update_check()
    assert started == []


# ── 2. default construction still starts the check ───────────────────────────

def test_default_construction_still_starts_the_check(app, started, monkeypatch):
    """The opt-out is additive: users still get the startup check."""
    monkeypatch.delenv(ENV, raising=False)
    win = mw.MainWindow()
    assert len(started) == 1
    assert isinstance(win._update_checker, UpdateChecker)
    assert started[0] is win._update_checker


def test_default_construction_still_shows_the_update_dialog(app, started, dialogs, monkeypatch):
    """...and an available update still raises the dialog for a real user."""
    monkeypatch.delenv(ENV, raising=False)
    win = mw.MainWindow()
    info = _info()
    win._on_update_available(info)
    assert dialogs == [info]


# ── 3. forcing update_available with the opt-out shows no dialog ─────────────

def test_forced_update_available_shows_no_dialog_when_opted_out(app, started, dialogs, monkeypatch):
    """The modal `dlg.exec()` an offscreen run cannot dismiss must never open."""
    monkeypatch.delenv(ENV, raising=False)
    win = mw.MainWindow(check_updates=False)
    win._on_update_available(_info())
    assert dialogs == []


def test_forced_update_available_shows_no_dialog_under_the_env_opt_out(app, started, dialogs):
    """Same, for the session-wide environment opt-out gui/conftest.py sets."""
    win = mw.MainWindow()
    win._on_update_available(_info())
    assert dialogs == []


# ── 4. Help → Check for Updates is unaffected ────────────────────────────────

def test_manual_check_still_runs_under_the_opt_out(app, started, monkeypatch):
    """The opt-out suppresses the *startup* check only.

    Help → Check for Updates is a deliberate user action; it must still fire, and
    must still be reachable when the window was built with `check_updates=False`.
    """
    win = mw.MainWindow(check_updates=False)
    win._manual_update_check()
    assert len(started) == 1
    assert started[0] is win._update_checker


def test_manual_check_still_shows_the_update_dialog_under_the_opt_out(app, started, dialogs):
    """And its own result handler still opens the dialog — nothing was weakened."""
    win = mw.MainWindow(check_updates=False)
    win._manual_update_check()
    info = _info()
    win._update_checker.update_available.emit(info)
    assert dialogs == [info]
