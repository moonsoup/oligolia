"""QThread workers for long-running bioinformatics operations."""

from __future__ import annotations
import sys
import os

# Ensure backend package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from PyQt6.QtCore import QThread, pyqtSignal
from typing import Any, Callable


class Worker(QThread):
    """Generic worker thread — runs any callable off the main thread."""
    result = pyqtSignal(object)
    error = pyqtSignal(str)
    progress = pyqtSignal(str)

    def __init__(self, fn: Callable, *args: Any, **kwargs: Any) -> None:
        super().__init__()
        self._fn = fn
        self._args = args
        self._kwargs = kwargs

    def run(self) -> None:
        try:
            result = self._fn(*self._args, **self._kwargs)
            self.result.emit(result)
        except Exception as e:
            self.error.emit(str(e))


class WorkerSlot:
    """Holds at most one live Worker, and refuses to start a second.

    #54: every panel kept its worker in a single attribute and overwrote it
    without checking whether the previous thread was still alive. Rebinding the
    attribute dropped the last reference to a *running* QThread, which Qt aborts
    on -- "QThread: Destroyed while thread is still running", exit 134, taking any
    unsaved sequences with it. Double-clicking Design Primers was enough.

    Six call sites had the identical bug, so the guard is here rather than
    repeated six times and fixed five of them.

    Refusing is deliberate rather than queueing: the panels are one-shot forms,
    and a queue would let a user stack up work they cannot see or cancel. A panel
    passes `on_busy` to say so in its own status line.
    """

    def __init__(self, on_busy: Callable[[], None] | None = None) -> None:
        self._worker: Worker | None = None
        self._on_busy = on_busy

    @property
    def worker(self) -> Worker | None:
        """The live worker, if any. Holding this is what keeps Qt from aborting."""
        return self._worker

    @property
    def busy(self) -> bool:
        return self._worker is not None and self._worker.isRunning()

    def start(
        self,
        fn: Callable,
        *args: Any,
        on_result: Callable[[Any], None] | None = None,
        on_error: Callable[[str], None] | None = None,
        on_progress: Callable[[str], None] | None = None,
        on_finished: Callable[[], None] | None = None,
        **kwargs: Any,
    ) -> Worker | None:
        """Run `fn` off the main thread, or return None if a job is already live.

        None means "refused, nothing started" -- the callable is not invoked.
        """
        if self.busy:
            if self._on_busy is not None:
                self._on_busy()
            return None

        worker = Worker(fn, *args, **kwargs)
        if on_result is not None:
            worker.result.connect(on_result)
        if on_error is not None:
            worker.error.connect(on_error)
        if on_progress is not None:
            worker.progress.connect(on_progress)
        if on_finished is not None:
            worker.finished.connect(on_finished)

        self._worker = worker
        worker.start()
        return worker

    def wait(self, msecs: int | None = None) -> bool:
        """Block until the current job finishes. For teardown and for tests."""
        if self._worker is None:
            return True
        return self._worker.wait() if msecs is None else self._worker.wait(msecs)


def worker_busy(owner: Any, *attrs: str) -> bool:
    """Is any of `owner`'s named worker attributes still running?

    The minimal guard for the existing panels (#54). Each of them already holds
    its worker in an attribute; the bug was never the holding, it was rebinding
    that attribute while the old QThread was still running, which dropped the
    last reference to a live thread and made Qt abort the process.

    So a call site needs one check before it constructs the next Worker, and
    nothing else about it has to change:

        if worker_busy(self, "_worker"):
            self._status.setText("Already working — wait for this one to finish.")
            return

    Several attributes may be named at once, for panels where two different
    buttons share one slot (pathways) or use separate ones (structure, search).

    `WorkerSlot` above is the tidier API and is what new panels should use; this
    exists so the fix to the existing six is one line each rather than a
    rewrite of six sets of signal wiring.
    """
    for attr in attrs:
        worker = getattr(owner, attr, None)
        if worker is not None and worker.isRunning():
            return True
    return False
