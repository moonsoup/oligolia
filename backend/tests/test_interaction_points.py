"""Tests for the putative-interaction-point heuristic (backend/services/interaction_points.py).

Uses a real, tiny reference structure (1CRN, crambin — a classic 46-residue test
protein) as a fixture rather than fabricated coordinates, so the SASA computation
is exercised against genuine atomic geometry.
"""

from pathlib import Path

from backend.services.interaction_points import compute_interaction_points

FIXTURE = Path(__file__).parent / "fixtures" / "1crn.pdb"


def _points() -> list[dict]:
    return compute_interaction_points(FIXTURE.read_text())


def test_returns_one_entry_per_standard_residue() -> None:
    points = _points()
    assert len(points) == 46  # crambin is a 46-residue protein


def test_exposed_acidic_residue_is_flagged() -> None:
    points = _points()
    asp43 = next(p for p in points if p["residue_index"] == 43)
    assert asp43["residue_name"] == "ASP"
    assert asp43["classification"] == "acidic"
    assert asp43["relative_sasa"] > 0.4
    assert asp43["is_putative_interaction_point"] is True


def test_buried_hydrophobic_residue_is_not_flagged() -> None:
    points = _points()
    phe13 = next(p for p in points if p["residue_index"] == 13)
    assert phe13["residue_name"] == "PHE"
    assert phe13["classification"] == "hydrophobic"
    assert phe13["relative_sasa"] < 0.1
    assert phe13["is_putative_interaction_point"] is False


def test_hydrophobic_residues_never_flagged_regardless_of_exposure() -> None:
    points = _points()
    for p in points:
        if p["classification"] == "hydrophobic":
            assert p["is_putative_interaction_point"] is False


def test_output_types_are_json_serializable() -> None:
    import json
    json.dumps(_points())  # raises on numpy scalar types
