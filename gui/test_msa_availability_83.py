"""#83 — the MSA tab says the aligner is missing before the user does the work.

Drives the real `AlignmentPanel` offscreen, the pattern CLAUDE.md prescribes:
build it in a `QApplication` under `QT_QPA_PLATFORM=offscreen`, exercise it
through its Qt API, assert on the resulting widget state and text.

Availability is controlled, never observed: `shutil.which` is patched to find
nothing, or a real (but inert) `muscle` executable is written into a tmp dir that
is put on PATH. A machine with MUSCLE actually installed runs exactly the same
assertions as one without.

The panel is built directly rather than through `MainWindow` — nothing here needs
the window, so nothing here needs #84's update-check opt-out either.
"""

from __future__ import annotations

import os
import shutil
import stat

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from gui.panels.alignment_panel import AlignmentPanel  # noqa: E402

#: Index of the MSA page in the panel's tab widget; the Pairwise tab is 0.
MSA_TAB = 1


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def no_aligners(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(shutil, "which", lambda *a, **k: None)


def _plant_muscle(directory, monkeypatch: pytest.MonkeyPatch) -> str:
    """Write an inert but executable `muscle` into `directory` and own PATH."""
    path = os.path.join(str(directory), "muscle")
    with open(path, "w") as fh:
        fh.write("#!/bin/sh\nexit 0\n")
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", str(directory))
    return path


def test_no_aligner_disables_the_msa_controls_and_names_the_install(
    app: QApplication, no_aligners: None, tmp_path
) -> None:
    """Acceptance 1 (panel half): disabled input and button, hint in the notice."""
    panel = AlignmentPanel()
    panel._tabs.setCurrentIndex(MSA_TAB)  # what the user does to reach the controls

    assert not panel._msa_input.isEnabled()
    assert not panel._btn_msa.isEnabled()

    notice = panel._msa_notice.text()
    assert panel._msa_notice.isVisibleTo(panel), "the notice must be showing"
    assert "brew install muscle" in notice, notice
    assert "drive5.com" in notice, "the hint's link is offered, not just its name"

    # The tab itself stays — #83 disables the controls, it does not remove the tab.
    assert panel._tabs.count() == 2
    assert "Multiple Sequence Alignment" in panel._tabs.tabText(MSA_TAB)

    panel.show()
    app.processEvents()
    panel.grab().save(str(tmp_path / "msa_unavailable.png"))


def test_fake_muscle_on_path_leaves_the_msa_controls_enabled(
    app: QApplication, tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Acceptance 2 (panel half): with an aligner found, the tab behaves as before."""
    _plant_muscle(tmp_path, monkeypatch)
    panel = AlignmentPanel()
    panel._tabs.setCurrentIndex(MSA_TAB)

    assert panel._msa_input.isEnabled()
    assert panel._btn_msa.isEnabled()
    assert not panel._msa_notice.isVisibleTo(panel), (
        f"nothing to install, so no notice: {panel._msa_notice.text()!r}"
    )


def test_installing_an_aligner_then_reshowing_the_tab_enables_the_controls(
    app: QApplication, tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Acceptance 3: unavailable → available → tab re-shown, same panel object."""
    monkeypatch.setattr(shutil, "which", lambda *a, **k: None)
    panel = AlignmentPanel()
    panel._tabs.setCurrentIndex(MSA_TAB)
    assert not panel._btn_msa.isEnabled()

    # The user installs MUSCLE while the app is running.
    monkeypatch.undo()
    path = _plant_muscle(tmp_path, monkeypatch)

    # Nothing re-checks on its own; the panel is not rebuilt.
    panel._tabs.setCurrentIndex(0)
    panel._tabs.setCurrentIndex(MSA_TAB)
    app.processEvents()

    assert panel._btn_msa.isEnabled(), "re-showing the tab must re-detect"
    assert panel._msa_input.isEnabled()
    assert not panel._msa_notice.isVisibleTo(panel)
    assert shutil.which("muscle") == path


def test_losing_the_aligner_disables_the_controls_again(
    app: QApplication, tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The re-check runs both ways, so a stale enabled button cannot linger."""
    _plant_muscle(tmp_path, monkeypatch)
    panel = AlignmentPanel()
    assert panel._btn_msa.isEnabled()

    monkeypatch.setattr(shutil, "which", lambda *a, **k: None)
    panel._tabs.setCurrentIndex(0)
    panel._tabs.setCurrentIndex(MSA_TAB)

    assert not panel._btn_msa.isEnabled()
    assert "brew install muscle" in panel._msa_notice.text()


def test_the_panel_and_the_router_read_the_same_detection(
    app: QApplication, no_aligners: None
) -> None:
    """The panel's notice is built from the router's hint table, not a copy of it."""
    from backend.routers.alignment import aligner_statuses

    panel = AlignmentPanel()
    muscle = next(a for a in aligner_statuses() if a.name == "muscle")
    assert muscle.available is False
    # Every word of advice the panel gives comes from the shared hint.
    assert muscle.hint.split(" (")[0] in panel._msa_notice.text()
    assert "brew install muscle" in muscle.hint


def test_pairwise_tab_is_untouched_when_no_aligner_is_installed(
    app: QApplication, no_aligners: None
) -> None:
    """MSA needs MUSCLE; pairwise never did, and must keep working without it."""
    panel = AlignmentPanel()
    panel._seq1.setPlainText("ATGGTGCACCTGACTCCTGAGGAGAAGTCT")
    panel._seq2.setPlainText("ATGGTGCACCTGACTCCTGAGGAGAAGTCA")

    panel._run_pairwise()

    shown = panel._pair_result.toPlainText()
    assert shown.startswith("Seq1"), shown
    assert "error" not in shown.lower(), shown
