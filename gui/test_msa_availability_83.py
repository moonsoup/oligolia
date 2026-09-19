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
from PyQt6.QtCore import Qt  # noqa: E402
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


# ── Round 2: the tab must run the aligner it says is available ───────────────
#
# Round 1 gated the controls on "any aligner available" but still ran the
# hard-coded default. On a ClustalW-only machine that meant enabled controls, no
# notice, and the missing-MUSCLE 503 after the user pasted sequences and pressed
# Run — the exact failure #83 exists to prevent, moved one aligner sideways. The
# selector is the fix, so these drive the selector.


def _plant(directory, monkeypatch: pytest.MonkeyPatch, *names: str) -> dict:
    """Put inert but executable copies of `names` on an otherwise empty PATH."""
    planted = {}
    for name in names:
        path = os.path.join(str(directory), name)
        with open(path, "w") as fh:
            fh.write("#!/bin/sh\nexit 0\n")
        os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        planted[name] = path
    monkeypatch.setenv("PATH", str(directory))
    return planted


def _record_runs(monkeypatch: pytest.MonkeyPatch) -> list:
    """Capture the `MSARequest`s the panel sends, without needing a real aligner.

    Patches the router module the panel calls through, so a panel that went back
    to calling its own private copy of MUSCLE would record nothing and fail.
    """
    from backend.routers import alignment

    calls: list = []

    def recording_multiple_align(req):
        calls.append(req)
        return alignment.MSAResult(
            aligned=[{"id": s["id"], "aligned_seq": s["seq"]} for s in req.sequences],
            consensus="N" * len(req.sequences[0]["seq"]),
            identity_matrix=[[100.0] * len(req.sequences)] * len(req.sequences),
        )

    monkeypatch.setattr(alignment, "multiple_align", recording_multiple_align)
    return calls


def _run_msa_and_wait(app: QApplication, panel: AlignmentPanel) -> None:
    """Click Run the way a user does, then let the worker thread finish."""
    panel._msa_input.setPlainText(">a\nATGGTGCACCTG\n>b\nATGGTGCATCTG")
    panel._btn_msa.click()
    assert panel._msa_worker is not None, "the click must have started a run"
    assert panel._msa_worker.wait(10_000), "the MSA worker did not finish"
    app.processEvents()


def test_clustalw_only_enables_the_controls_and_runs_clustalw(
    app: QApplication, tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The verifier's round-1 counterexample: ClustalW installed, MUSCLE not.

    The controls are enabled — ClustalW can align these sequences — and the Run
    that follows must ask for ClustalW. Running the hard-coded MUSCLE here is
    what produced a 503 after the work, which is the bug #83 is about.
    """
    _plant(tmp_path, monkeypatch, "clustalw")
    calls = _record_runs(monkeypatch)

    panel = AlignmentPanel()
    panel._tabs.setCurrentIndex(MSA_TAB)

    assert panel._msa_input.isEnabled()
    assert panel._btn_msa.isEnabled()
    assert not panel._msa_notice.isVisibleTo(panel)
    assert panel._msa_algorithm.currentData() == "clustalw"
    assert "ClustalW" in panel._btn_msa.text(), panel._btn_msa.text()

    _run_msa_and_wait(app, panel)

    assert [c.algorithm for c in calls] == ["clustalw"]
    assert "Error" not in panel._msa_status.text(), panel._msa_status.text()


def test_muscle_is_preferred_and_run_when_both_are_installed(
    app: QApplication, tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With everything installed the tab behaves as it always did: MUSCLE."""
    _plant(tmp_path, monkeypatch, "muscle", "clustalw")
    calls = _record_runs(monkeypatch)

    panel = AlignmentPanel()
    panel._tabs.setCurrentIndex(MSA_TAB)
    assert panel._msa_algorithm.currentData() == "muscle"

    _run_msa_and_wait(app, panel)
    assert [c.algorithm for c in calls] == ["muscle"]


def test_choosing_the_other_aligner_runs_that_one(
    app: QApplication, tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The selector is not decoration — what it shows is what gets run."""
    _plant(tmp_path, monkeypatch, "muscle", "clustalw")
    calls = _record_runs(monkeypatch)

    panel = AlignmentPanel()
    panel._tabs.setCurrentIndex(MSA_TAB)
    panel._msa_algorithm.setCurrentIndex(panel._msa_algorithm.findData("clustalw"))

    _run_msa_and_wait(app, panel)
    assert [c.algorithm for c in calls] == ["clustalw"]


def test_a_missing_aligner_is_listed_but_cannot_be_chosen(
    app: QApplication, tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """MUSCLE stays visible while missing — greyed, hinted, not selectable."""
    _plant(tmp_path, monkeypatch, "clustalw")

    panel = AlignmentPanel()
    panel._tabs.setCurrentIndex(MSA_TAB)

    names = [panel._msa_algorithm.itemData(i) for i in range(panel._msa_algorithm.count())]
    assert names == ["muscle", "clustalw"], names

    model = panel._msa_algorithm.model()
    assert not model.item(names.index("muscle")).isEnabled()
    assert model.item(names.index("clustalw")).isEnabled()
    assert "not installed" in panel._msa_algorithm.itemText(names.index("muscle")).lower()
    assert "brew install muscle" in panel._msa_algorithm.itemData(
        names.index("muscle"), Qt.ItemDataRole.ToolTipRole
    )


def test_no_aligner_leaves_nothing_selected_and_nothing_runnable(
    app: QApplication, no_aligners: None
) -> None:
    """With nothing installed the selector is disabled and picks nothing."""
    panel = AlignmentPanel()
    panel._tabs.setCurrentIndex(MSA_TAB)

    assert not panel._msa_algorithm.isEnabled()
    assert panel._msa_algorithm.currentData() is None
    assert panel._btn_msa.text() == "Run MSA (requires MUSCLE)"


def test_the_selector_offers_exactly_what_the_router_reports(
    app: QApplication, tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One decision: the panel's default is the router's `preferred`."""
    from backend.routers.alignment import preferred_algorithm

    _plant(tmp_path, monkeypatch, "clustalw")
    panel = AlignmentPanel()
    assert panel._msa_algorithm.currentData() == preferred_algorithm() == "clustalw"


def test_losing_the_chosen_aligner_falls_back_on_re_show(
    app: QApplication, tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A pick that stops being installed is replaced, not kept and then 503'd."""
    _plant(tmp_path, monkeypatch, "muscle", "clustalw")
    panel = AlignmentPanel()
    panel._tabs.setCurrentIndex(MSA_TAB)
    panel._msa_algorithm.setCurrentIndex(panel._msa_algorithm.findData("clustalw"))

    # ClustalW is uninstalled while the app is open; MUSCLE stays.
    os.unlink(os.path.join(str(tmp_path), "clustalw"))
    panel._tabs.setCurrentIndex(0)
    panel._tabs.setCurrentIndex(MSA_TAB)

    assert panel._msa_algorithm.currentData() == "muscle"
    assert panel._btn_msa.isEnabled()
