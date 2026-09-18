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
    assert shown.startswith("Seq1"), shown
    assert len(shown) > 50, shown
    # Identity/score live in the stat labels above the pane, not duplicated in it (#77).
    assert "—" not in panel._stat_identity.text(), panel._stat_identity.text()
    assert "%" in panel._stat_identity.text()
    assert "—" not in panel._stat_score.text()


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


# --- #77: the alignment view must line up ---

from gui.panels.alignment_panel import BLOCK_WIDTH, format_alignment_blocks  # noqa: E402


def _rows(text: str) -> list[list[str]]:
    """Group the rendered output back into (seq1, match, seq2) triples."""
    lines = text.split("\n")
    blocks = []
    for i in range(0, len(lines), 4):
        triple = lines[i:i + 3]
        if len(triple) == 3 and triple[0].startswith("Seq1"):
            blocks.append(triple)
    return blocks


def test_every_block_row_has_the_same_width() -> None:
    """#77: the three rows wrapped at different columns, so the bars drifted."""
    a1 = "ACGTACGTAC" * 12
    a2 = "ACGTAGGTAC" * 12
    out = format_alignment_blocks(a1, a2)
    blocks = _rows(out)
    assert blocks, out
    for seq1, match, seq2 in blocks:
        assert len(seq1) == len(seq2), (seq1, seq2)
        # Strictly equal, not merely "not longer" — Codex noted the looser form let
        # a short match row through while the commit claimed equal widths.
        assert len(match) == len(seq1), (len(match), len(seq1), match, seq1)


def test_every_bar_sits_under_two_matching_bases() -> None:
    """The property the old rendering violated, stated directly."""
    a1 = "ACGT-ACGTACGTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
    a2 = "ACGTAACGTAGGTTTACGTACG--CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
    out = format_alignment_blocks(a1, a2)

    for seq1_row, match_row, seq2_row in _rows(out):
        # Strip the identical label+coordinate prefix from each row.
        prefix = len(seq1_row) - len(seq1_row.lstrip())  # not reliable; use fixed parse
        parts1 = seq1_row.split()
        parts2 = seq2_row.split()
        chunk1, chunk2 = parts1[2], parts2[2]
        offset = seq1_row.index(chunk1)
        match = match_row[offset:offset + len(chunk1)]
        # No `or match.strip() == ""` escape: Codex pointed out that let a blank
        # match row satisfy the test, which is exactly the vacuous case.
        assert len(match) == len(chunk1), (len(match), len(chunk1))

        for i, ch in enumerate(match):
            if ch == "|":
                assert chunk1[i] == chunk2[i] != "-", (i, chunk1[i], chunk2[i])
            elif ch == ".":
                assert chunk1[i] != chunk2[i]
                assert chunk1[i] != "-" and chunk2[i] != "-"
        del prefix


def test_coordinates_count_bases_not_columns() -> None:
    """A gap must not advance the coordinate, or it cannot locate a mismatch."""
    a1 = "ACGT" + "-" * 10 + "ACGT" * 20
    a2 = "ACGT" + "ACGTACGTAC" + "ACGT" * 20
    out = format_alignment_blocks(a1, a2)
    first = _rows(out)[0]
    # Seq1's end coordinate for the first block excludes the 10 gap characters.
    end1 = int(first[0].split()[-1])
    end2 = int(first[2].split()[-1])
    assert end1 == BLOCK_WIDTH - 10, (end1, BLOCK_WIDTH)
    assert end2 == BLOCK_WIDTH, end2


def test_mismatched_lengths_are_refused() -> None:
    """Don't render something that cannot line up."""
    with _pytest.raises(ValueError):
        format_alignment_blocks("ACGT", "ACG")


def test_blocks_are_no_wider_than_the_block_width() -> None:
    a1 = "ACGT" * 50
    a2 = "ACGA" * 50
    for seq1_row, _m, _s2 in _rows(format_alignment_blocks(a1, a2)):
        assert len(seq1_row.split()[2]) <= BLOCK_WIDTH


def test_the_panel_uses_the_block_formatter_and_does_not_wrap() -> None:
    """Codex: a substring search for the name also matches its own definition.

    So resolve it in the AST and require the call to happen inside _run_pairwise
    — otherwise the pin passes even if that method goes back to rendering the
    alignment by hand.
    """
    import ast
    import inspect
    import textwrap

    from gui.panels.alignment_panel import AlignmentPanel

    src = textwrap.dedent(inspect.getsource(AlignmentPanel._run_pairwise))
    fn = ast.parse(src).body[0]
    called = {
        n.func.id for n in ast.walk(fn)
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
    }
    assert "format_alignment_blocks" in called, (
        "_run_pairwise must render through format_alignment_blocks (#77)"
    )

    from pathlib import Path

    source = Path(__file__).with_name("alignment_panel.py").read_text()
    assert "LineWrapMode.NoWrap" in source, "the result pane must not word-wrap (#77)"
