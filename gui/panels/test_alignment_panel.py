"""Drive the real AlignmentPanel widget offscreen.

The pattern CLAUDE.md prescribes for GUI changes: build the widget in a real
QApplication under QT_QPA_PLATFORM=offscreen, exercise it through its Qt API, and
assert on the resulting state rather than on internals.

#56 showed up in this panel as a blank "Alignment error: " for ordinary input,
because `list(aligner.align(...))` tried to materialise more co-optimal alignments
than fit in an int64. Nothing in the backend suite could see that — the panel does
its own aligning inline rather than calling the router.
"""

from __future__ import annotations

import os
import random

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from gui.panels.alignment_panel import AlignmentPanel  # noqa: E402


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


def _unrelated(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_pairwise_on_unrelated_500nt_shows_a_result_not_an_error(app: QApplication) -> None:
    """#56: the user-visible symptom was `Alignment error: ` with no message."""
    panel = AlignmentPanel()
    panel._seq1.setPlainText(_unrelated(500, 1))
    panel._seq2.setPlainText(_unrelated(500, 2))

    panel._run_pairwise()

    shown = panel._pair_result.toPlainText()
    assert "error" not in shown.lower(), shown
    assert "No alignment found" not in shown, shown
    # The panel reports identity/score; both must be present and parseable.
    assert "Identity" in shown or "identity" in shown, shown
    assert len(shown) > 50, shown


# Deliberately NOT tested here:
#   - 2 kb x 2 kb. The panel renders the full gapped alignment, which takes ~2 minutes
#     at that size; `backend/tests/test_alignment.py` covers 2 kb at the router level in
#     well under a second. A two-minute GUI test is one nobody will run.
#   - empty input. The panel guards it correctly (`if not s1 or not s2`) but reports it
#     with a modal `QMessageBox.warning`, which blocks forever with no one to dismiss
#     it. Driving that path headlessly needs the dialog patched out; the guard itself is
#     plain enough to read.


def test_screenshot_of_a_completed_alignment(app: QApplication, tmp_path) -> None:
    """Render the panel after a real alignment, so a human can look at it."""
    panel = AlignmentPanel()
    panel.resize(1000, 700)
    panel._seq1.setPlainText(_unrelated(300, 5))
    panel._seq2.setPlainText(_unrelated(300, 6))
    panel._run_pairwise()

    out = tmp_path / "alignment_panel.png"
    assert panel.grab().save(str(out)), "grab().save() failed"
    assert out.stat().st_size > 1000, "screenshot suspiciously small"


# --- #58: the panel must not fabricate an alignment either ---

import shutil  # noqa: E402

import pytest as _pytest  # noqa: E402


def test_msa_without_an_aligner_raises_instead_of_padding(app: QApplication) -> None:
    """#58: the panel had its own copy of the right-padding fallback.

    Worker.error is already wired to the status line, so raising is enough to tell
    the user the truth. What must never happen is a dict of padded input coming
    back and `_on_msa_done` reporting "Aligned 2 sequences."
    """
    if shutil.which("muscle"):
        _pytest.skip("muscle is installed here; this test is about its absence")

    panel = AlignmentPanel()
    seqs = [
        {"id": "a", "seq": "ATGGTGCACCTGACTCCTGAGGAGAAGTCT"},
        {"id": "b", "seq": "GATGGTGCACCTGACTCCTGAGGAGAAGTCT"},
    ]

    with _pytest.raises(RuntimeError) as err:
        panel._do_msa(seqs)

    msg = str(err.value).lower()
    assert "muscle" in msg
    assert "not installed" in msg


def test_the_panel_no_longer_carries_its_own_aligner_copy() -> None:
    """#58 also asked for the duplication to go: one implementation, one fix."""
    from pathlib import Path

    source = Path(__file__).with_name("alignment_panel.py").read_text()
    assert "muscle" not in source.replace('algorithm="muscle"', ""), (
        "the panel should delegate to backend.routers.alignment, not shell out itself (#58)"
    )
    assert "multiple_align" in source
