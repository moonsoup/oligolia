"""#64.4: the CRISPR panel must not quietly shrink the search space.

`set_target` did `seq[:2000]` with the comment "limit to 2 kbp for display", but
`_run_design` reads the box — so the truncation changed the *results*, not the
display. A 5 kb insert was designed against its first 2 kb and nothing said so.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from gui.panels.crispr_panel import CRISPRPanel  # noqa: E402


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


def _seq(n: int) -> str:
    return ("ACGTGGCATCACGATGGCCTGTGGG" * (n // 25 + 1))[:n]


def test_a_long_target_is_loaded_whole(app: QApplication) -> None:
    panel = CRISPRPanel()
    seq = _seq(5000)

    panel.set_target(seq)

    loaded = panel._target_input.toPlainText()
    assert len(loaded) == 5000, f"target was truncated to {len(loaded)} nt (#64.4)"
    assert loaded == seq


def test_a_long_target_says_so(app: QApplication) -> None:
    panel = CRISPRPanel()
    panel.set_target(_seq(5000))
    status = panel._status.text()
    assert "5,000" in status, status
    assert "whole sequence" in status.lower(), status


def test_a_short_target_gets_no_notice(app: QApplication) -> None:
    panel = CRISPRPanel()
    panel.set_target(_seq(300))
    assert panel._target_input.toPlainText() == _seq(300)
    assert panel._status.text() == ""


def test_set_target_contains_no_slice_at_all(app: QApplication) -> None:
    """Pin it in the AST, so the docstring that explains the history cannot trip it.

    Any slice of the incoming sequence inside set_target is the defect, whatever
    bound it uses — 2000 was just the number that happened to be there.
    """
    import ast
    import inspect
    import textwrap

    src = textwrap.dedent(inspect.getsource(CRISPRPanel.set_target))
    fn = ast.parse(src).body[0]
    slices = [n for n in ast.walk(fn) if isinstance(n, ast.Subscript) and isinstance(n.slice, ast.Slice)]
    assert not slices, "set_target must not truncate the target (#64.4)"


# --- #64 interaction: the panel must not pin guide_length ---

def test_the_panel_does_not_pass_guide_length(app: QApplication) -> None:
    """Caught while fixing #54, and it is a real #64 regression.

    The router honours `guide_length` only when a caller explicitly sets it, and
    otherwise uses each nuclease's canonical length. The panel was passing
    `guide_length=20` unconditionally, so after #64 every Cas12a guide designed
    from the GUI would have been 20 nt instead of 23 — a silent shortening that
    no backend test could see, because the backend was being asked for 20.
    """
    import ast
    import inspect
    import textwrap

    src = textwrap.dedent(inspect.getsource(CRISPRPanel._run_design))
    fn = ast.parse(src).body[0]
    for node in ast.walk(fn):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) \
                and node.func.id == "CRISPRDesignRequest":
            names = {kw.arg for kw in node.keywords}
            assert "guide_length" not in names, (
                "the panel must not pin guide_length; the router picks the "
                "nuclease's canonical length (#64)"
            )
            break
    else:
        raise AssertionError("CRISPRDesignRequest(...) not found in _run_design")


def test_cas12a_from_the_panels_request_gets_23nt_guides(app: QApplication) -> None:
    """End to end through the router, with the panel's own request shape."""
    from backend.models.crispr import CasType, CRISPRDesignRequest
    from backend.routers.crispr import design_guides

    target = "TTTA" + "ACGTGGCATC" * 10
    req = CRISPRDesignRequest(
        target_sequence=target, cas_type=CasType.CAS12A, max_guides=10,
        check_off_targets=False,
    )
    for g in design_guides(req).guides:
        assert len(g.sequence) == 23, g.sequence
