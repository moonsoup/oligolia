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
