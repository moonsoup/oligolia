"""Tests for backend/routers/structure.py — network calls are always mocked."""

from pathlib import Path

from backend.services.pdb import PDBClient
from backend.services.alphafold_db import AlphaFoldDBClient, AlphaFoldPrediction
from backend.services.esmfold import ESMFoldClient, MAX_LENGTH

FIXTURE = Path(__file__).parent / "fixtures" / "1crn.pdb"


def test_experimental_hit_is_preferred(client, monkeypatch) -> None:
    monkeypatch.setattr(PDBClient, "search_by_uniprot", lambda self, uid, max_results=10: ["1CRN"])
    monkeypatch.setattr(PDBClient, "download_pdb", lambda self, pdb_id: FIXTURE.read_text())

    resp = client.post("/structure/predict", json={"sequence": "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN", "uniprot_id": "P01542"})
    assert resp.status_code == 200
    data = resp.json()
    assert data["source"] == "experimental_pdb"
    assert data["pdb_id"] == "1CRN"


def test_falls_back_to_alphafold_db_when_no_experimental_hit(client, monkeypatch) -> None:
    monkeypatch.setattr(PDBClient, "search_by_uniprot", lambda self, uid, max_results=10: [])
    monkeypatch.setattr(
        AlphaFoldDBClient, "get_prediction",
        lambda self, uid: AlphaFoldPrediction(
            uniprot_accession=uid, pdb_text=FIXTURE.read_text(), global_metric_value=88.0,
        ),
    )

    resp = client.post("/structure/predict", json={"sequence": "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN", "uniprot_id": "P01542"})
    assert resp.status_code == 200
    data = resp.json()
    assert data["source"] == "predicted_alphafold_db"
    assert "88.0" in data["confidence_note"]


def test_falls_back_to_esmfold_when_no_uniprot_id(client, monkeypatch) -> None:
    monkeypatch.setattr(ESMFoldClient, "predict", lambda self, seq: FIXTURE.read_text())

    resp = client.post("/structure/predict", json={"sequence": "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"})
    assert resp.status_code == 200
    data = resp.json()
    assert data["source"] == "predicted_esmfold"


def test_over_length_sequence_gives_400_not_a_crash(client) -> None:
    too_long = "A" * (MAX_LENGTH + 1)
    resp = client.post("/structure/predict", json={"sequence": too_long})
    assert resp.status_code == 400
    assert str(MAX_LENGTH) in resp.json()["detail"]


def test_interaction_points_endpoint_shape(client) -> None:
    resp = client.post("/structure/interaction_points", json={"pdb_text": FIXTURE.read_text()})
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["points"]) == 46
    assert "heuristic" in data["method"].lower()


def test_interaction_points_endpoint_rejects_garbage_input(client) -> None:
    resp = client.post("/structure/interaction_points", json={"pdb_text": "not a pdb file"})
    assert resp.status_code == 400


# --- #66.3: say which PDB entry was chosen, and on what basis ---

def test_the_chosen_entry_and_its_alternatives_are_reported(monkeypatch) -> None:
    """`ids[0]` from an unranked search, with a note that mentioned neither.

    Full sequence-identity ranking needs a live fetch per candidate and is
    follow-up work; what this fixes is that the user could not tell an arbitrary
    pick from a good one.
    """
    from backend.models.structure import StructureRequest
    from backend.routers import structure as mod

    monkeypatch.setattr(mod._pdb, "search_by_uniprot", lambda *a, **k: ["1ABC", "2DEF", "3GHI"])
    monkeypatch.setattr(mod._pdb, "download_pdb", lambda pdb_id: f"HEADER {pdb_id}\nEND\n")

    result = mod.get_or_predict_structure(
        StructureRequest(sequence="MVHLTPEEK", uniprot_id="P68871")
    )

    assert result.pdb_id == "1ABC"
    assert result.pdb_candidates == ["1ABC", "2DEF", "3GHI"]
    note = result.confidence_note
    assert "2DEF" in note, note
    assert "3 entries" in note, note
    # It must not imply an identity or organism check that did not happen.
    assert "Not ranked by sequence identity" in note, note
    assert "not" in note.lower() and "organism" in note.lower(), note


def test_a_single_hit_does_not_claim_alternatives(monkeypatch) -> None:
    from backend.models.structure import StructureRequest
    from backend.routers import structure as mod

    monkeypatch.setattr(mod._pdb, "search_by_uniprot", lambda *a, **k: ["1ABC"])
    monkeypatch.setattr(mod._pdb, "download_pdb", lambda pdb_id: "HEADER\nEND\n")

    result = mod.get_or_predict_structure(
        StructureRequest(sequence="MVHLTPEEK", uniprot_id="P68871")
    )
    assert result.pdb_candidates == ["1ABC"]
    assert "out of" not in result.confidence_note
