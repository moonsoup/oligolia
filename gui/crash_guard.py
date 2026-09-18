"""Keep an unhandled exception from killing the app.

PyQt6 aborts the process (exit -6) when an exception escapes a slot, so without
this every missed edge case is a crash with nothing on screen and nothing on disk.
Confirmed paths at the time of writing (#55): saving FASTA or GenBank to an
unwritable location, and entering `[]` as workflow step parameters, where only
`json.JSONDecodeError` was caught and pydantic's `ValidationError` escaped.

This is the general guard. Specific handlers still belong at their call sites — a
dialog saying "PermissionError" is much worse than one saying "that folder is not
writable" — but a missing handler should now cost the user a dialog, not their
session.

Deliberately: `KeyboardInterrupt` is delegated to the default hook, because Ctrl-C
must still stop the program; swallowing it would be worse than crashing. And every
step is individually guarded, because this code runs while the interpreter is
already unwinding — a crash guard that raises is just a second crash.
"""

from __future__ import annotations

import sys
import tempfile
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

#: Same convention as gui/updater.py's oligolia_update.log.
LOG_PATH = Path(tempfile.gettempdir()) / "oligolia_crash.log"

TITLE = "Oligolia hit an unexpected error"


def _default_show(title: str, body: str) -> None:
    """Show the dialog. Imported lazily so this module is usable without Qt up."""
    from PyQt6.QtWidgets import QApplication, QMessageBox

    if QApplication.instance() is None:
        # No event loop to attach a dialog to; the log and stderr are all we have.
        print(body, file=sys.stderr)
        return
    QMessageBox.critical(None, title, body)


def handle_exception(
    exc_type: type[BaseException],
    exc: BaseException,
    tb: object,
    *,
    show: Callable[[str, str], None] | None = None,
    log_path: Path | None = None,
    fallback: Callable[..., None] | None = None,
) -> None:
    """Log the traceback, tell the user, and return so the app keeps running."""
    if issubclass(exc_type, KeyboardInterrupt):
        (fallback or sys.__excepthook__)(exc_type, exc, tb)  # type: ignore[arg-type]
        return

    path = log_path or LOG_PATH
    detail = "".join(traceback.format_exception(exc_type, exc, tb))  # type: ignore[arg-type]

    # Always get it to stderr, which needs nothing to work.
    try:
        print(detail, file=sys.stderr)
    except Exception:
        pass

    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8") as fh:
            fh.write(f"\n=== {datetime.now(timezone.utc).isoformat()} ===\n")
            fh.write(detail)
    except Exception:
        pass  # an unwritable log must not become the crash

    body = (
        f"{exc_type.__name__}: {exc}\n\n"
        "The app is still running, but this action did not finish. "
        "If it keeps happening, the full details are in:\n"
        f"{path}"
    )
    try:
        (show or _default_show)(TITLE, body)
    except Exception:
        pass  # no display, or Qt already tearing down


def install(show: Callable[[str, str], None] | None = None) -> Callable[[], None]:
    """Install the hook. Returns a callable that restores the previous one."""
    previous = sys.excepthook

    def hook(exc_type, exc, tb):  # noqa: ANN001, ANN202
        handle_exception(exc_type, exc, tb, show=show, fallback=previous)

    sys.excepthook = hook
    return lambda: setattr(sys, "excepthook", previous)
