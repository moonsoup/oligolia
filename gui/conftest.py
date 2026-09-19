"""Session-wide guards for the offscreen GUI suite.

#84: `MainWindow.__init__` starts a live `UpdateChecker` — an outbound request to
the GitHub releases API from a QThread — and, if the code under test is older than
the latest published release, `update_available` fires and opens a *modal* dialog
that nothing can dismiss under `QT_QPA_PLATFORM=offscreen`. That is exactly the
situation when a verifier re-runs an older commit after the next release tag, so
the suite would hang rather than fail.

`MainWindow(check_updates=False)` is the explicit opt-out and the existing tests
use it. This file is the belt to that pair of braces: setting the variable at
import time, before any test module is collected, means a GUI test written later
makes no request even if its author never heard of the flag. Tests that need the
default behaviour back (`gui/test_update_check_84.py` pins it) delete the variable
with `monkeypatch.delenv`, which restores it at teardown.
"""

from __future__ import annotations

import os

# Before importing anything that pulls in PyQt6, as every test module here does.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from gui.main_window import NO_UPDATE_CHECK_ENV  # noqa: E402

os.environ[NO_UPDATE_CHECK_ENV] = "1"
