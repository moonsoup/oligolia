"""#66.1: SASA was computed across every model at once, so they occluded each other.

An NMR ensemble holds many models of the SAME molecule at slightly different
conformations, all in one coordinate frame. Running Shrake-Rupley over the whole
structure lets them bury each other's surface: duplicating 1CRN's single model
dropped flagged interaction points 15 -> 8 and mean relative SASA 0.376 -> 0.135.

HETATM waters occlude by the same path.

Uses the committed 1CRN fixture, so this runs offline.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from backend.services.interaction_points import compute_interaction_points

FIXTURE = Path(__file__).parent / "fixtures" / "1crn.pdb"

pytestmark = pytest.mark.skipif(not FIXTURE.is_file(), reason="1crn.pdb fixture missing")


def _single() -> str:
    return FIXTURE.read_text()


def _duplicated_model() -> str:
    """The same molecule as two models, the shape of an NMR ensemble."""
    lines = [
        ln for ln in _single().splitlines()
        if not ln.startswith(("MODEL", "ENDMDL", "END", "MASTER", "CONECT"))
    ]
    out = ["MODEL        1"] + lines + ["ENDMDL", "MODEL        2"] + lines + ["ENDMDL", "END"]
    return "\n".join(out) + "\n"


def _with_waters() -> str:
    """The same molecule with waters packed right against the surface."""
    lines = _single().splitlines()
    waters = [
        f"HETATM{9000 + i:5d}  O   HOH A{500 + i:4d}    "
        f"{10.0 + i * 0.5:8.3f}{10.0:8.3f}{10.0:8.3f}  1.00 20.00           O"
        for i in range(40)
    ]
    # Immediately after the last ATOM record. Placing them after TER/MASTER
    # meant the parser never took them, which made this test vacuous.
    last_atom = max(i for i, ln in enumerate(lines) if ln.startswith("ATOM"))
    return "\n".join(lines[: last_atom + 1] + waters + lines[last_atom + 1 :]) + "\n"


def _flagged(points: list[dict]) -> int:
    return sum(1 for p in points if p["is_putative_interaction_point"])


def _mean_sasa(points: list[dict]) -> float:
    return sum(p["relative_sasa"] for p in points) / len(points)


def test_a_duplicated_model_does_not_change_the_answer(app=None) -> None:
    """The headline number from the issue: 15 flagged points became 8."""
    one = compute_interaction_points(_single())
    two = compute_interaction_points(_duplicated_model())

    assert _flagged(one) == _flagged(two), (
        f"{_flagged(one)} interaction points with one model, {_flagged(two)} with "
        "the same molecule duplicated — the models are occluding each other"
    )
    assert abs(_mean_sasa(one) - _mean_sasa(two)) < 0.01, (
        _mean_sasa(one), _mean_sasa(two)
    )


def test_only_one_model_is_reported() -> None:
    """Guard on behaviour that was already correct.

    The reporting loop always took the first model (there was a `break` for it),
    which is why the output LOOKED right: the residue list was from model 1 and
    only its SASA numbers had been occluded by the others. Worth pinning anyway,
    since the model selection moved when the SASA fix landed.
    """
    one = compute_interaction_points(_single())
    two = compute_interaction_points(_duplicated_model())
    assert len(two) == len(one), (len(one), len(two))

    keys = [(p["chain"], p["residue_index"]) for p in two]
    assert len(keys) == len(set(keys)), "a residue is reported more than once"


def test_waters_do_not_bury_the_surface() -> None:
    """HETATM waters occlude by the same path as extra models."""
    bare = compute_interaction_points(_single())
    wet = compute_interaction_points(_with_waters())

    assert _flagged(bare) == _flagged(wet), (
        f"{_flagged(bare)} interaction points without waters, {_flagged(wet)} with "
        "them — waters are being counted as occluding atoms"
    )
    assert len(wet) == len(bare), "waters leaked into the residue list"


def test_the_single_model_answer_is_unchanged() -> None:
    """Regression guard: the fix must not alter the ordinary single-model result.

    1CRN is a 46-residue crambin structure with one model and no waters, so it
    is the case where nothing should move.
    """
    points = compute_interaction_points(_single())
    assert len(points) == 46, len(points)
    assert _flagged(points) > 0
    assert 0.0 < _mean_sasa(points) < 1.0
