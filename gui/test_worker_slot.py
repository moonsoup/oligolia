"""#54: starting a second job while one runs aborted the process.

"QThread: Destroyed while thread is still running", exit 134, taking unsaved
sequences with it. Every panel kept its worker in a single attribute and
overwrote it without checking whether the previous thread was still alive, so the
old QThread was garbage-collected mid-run.

Six call sites had the same bug, so the guard lives in one place.
"""

from __future__ import annotations

import os
import time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from gui.workers import WorkerSlot  # noqa: E402


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


def _slow(seconds: float = 0.4) -> str:
    time.sleep(seconds)
    return "done"


def test_a_fresh_slot_is_not_busy(app) -> None:
    assert WorkerSlot().busy is False


def test_starting_a_job_makes_the_slot_busy(app) -> None:
    slot = WorkerSlot()
    worker = slot.start(_slow, 0.3)
    try:
        assert worker is not None
        assert slot.busy is True
    finally:
        slot.wait()


def test_a_second_start_while_busy_is_refused(app) -> None:
    """The crash: the second start replaced the attribute and freed a live QThread."""
    slot = WorkerSlot()
    first = slot.start(_slow, 0.4)
    try:
        second = slot.start(_slow, 0.4)
        assert second is None, "a second worker was started while the first was live"
        # And the first must still be the one we hold.
        assert slot.worker is first
        assert first.isRunning()
    finally:
        slot.wait()


def test_the_slot_is_reusable_once_the_job_finishes(app) -> None:
    slot = WorkerSlot()
    slot.start(_slow, 0.05)
    slot.wait()
    assert slot.busy is False
    assert slot.start(_slow, 0.05) is not None
    slot.wait()


def test_the_worker_reference_is_held_for_the_whole_run(app) -> None:
    """The mechanism of the abort: nothing kept the QThread alive."""
    slot = WorkerSlot()
    worker = slot.start(_slow, 0.3)
    try:
        assert slot.worker is worker, "the slot must hold the worker while it runs"
        assert worker.isRunning()
    finally:
        slot.wait()


def test_a_refused_start_does_not_run_the_callable(app) -> None:
    calls: list[int] = []

    def record() -> None:
        calls.append(1)
        time.sleep(0.3)

    slot = WorkerSlot()
    slot.start(record)
    try:
        slot.start(record)
        slot.start(record)
    finally:
        slot.wait()

    assert calls == [1], f"the callable ran {len(calls)} times; only the first start should run"


def test_on_busy_is_called_so_a_panel_can_tell_the_user(app) -> None:
    told: list[str] = []
    slot = WorkerSlot(on_busy=lambda: told.append("busy"))
    slot.start(_slow, 0.3)
    try:
        slot.start(_slow, 0.3)
        assert told == ["busy"]
    finally:
        slot.wait()


def test_wait_returns_promptly_on_an_idle_slot(app) -> None:
    started = time.monotonic()
    WorkerSlot().wait()
    assert time.monotonic() - started < 1.0


# --- worker_busy: the minimal guard used by the existing panels ---

from gui.workers import Worker, worker_busy  # noqa: E402


class _Owner:
    pass


def test_worker_busy_is_false_when_the_attribute_is_absent(app) -> None:
    assert worker_busy(_Owner(), "_worker") is False


def test_worker_busy_is_false_when_the_attribute_is_none(app) -> None:
    owner = _Owner()
    owner._worker = None
    assert worker_busy(owner, "_worker") is False


def test_worker_busy_is_true_while_a_thread_runs(app) -> None:
    owner = _Owner()
    owner._worker = Worker(_slow, 0.3)
    owner._worker.start()
    try:
        assert worker_busy(owner, "_worker") is True
    finally:
        owner._worker.wait()
    assert worker_busy(owner, "_worker") is False


def test_worker_busy_checks_every_named_attribute(app) -> None:
    """Panels with two buttons keep two attributes; either being live counts."""
    owner = _Owner()
    owner._worker = None
    owner._points_worker = Worker(_slow, 0.3)
    owner._points_worker.start()
    try:
        assert worker_busy(owner, "_worker", "_points_worker") is True
        assert worker_busy(owner, "_worker") is False
    finally:
        owner._points_worker.wait()


# --- the live repro from the issue ---

_DOUBLE_CLICK = '''
import os, sys
os.environ["QT_QPA_PLATFORM"] = "offscreen"
sys.path.insert(0, {root!r})

from PyQt6.QtCore import QTimer
from PyQt6.QtWidgets import QApplication

app = QApplication([])
from gui.panels.crispr_panel import CRISPRPanel

BIG = ("ATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACGTGGCATCACGTTT"
       "CGGCCATGGACGTAGGCATCACGTGGCATCACGTGGCCGGCGGCGGCGCTGCTGCAAGCTT") * 160

panel = CRISPRPanel()
panel.set_target(BIG)
panel._run_design()
panel._run_design()
panel._run_design()

QTimer.singleShot(1500, app.quit)
sys.exit(app.exec())
'''


def test_clicking_design_three_times_does_not_abort_the_process(tmp_path) -> None:
    """#54's repro, run for real.

    Measured on this machine: without the guard the process ends with exit 134 and
    "QThread: Destroyed while thread is still running"; with it, exit 0. Verified
    both ways by temporarily removing the guard, so this test has teeth rather
    than asserting that a passing thing passes.
    """
    import subprocess
    import sys as _sys
    from pathlib import Path

    root = str(Path(__file__).resolve().parents[1])
    script = tmp_path / "double_click.py"
    script.write_text(_DOUBLE_CLICK.format(root=root))

    done = subprocess.run(
        [_sys.executable, str(script)], capture_output=True, text=True, timeout=180
    )
    assert done.returncode == 0, (
        f"exit {done.returncode} — 134/-6 means the QThread abort is back (#54)\n"
        f"{done.stderr[-2000:]}"
    )
    assert "Destroyed while thread is still running" not in done.stderr, done.stderr[-2000:]
