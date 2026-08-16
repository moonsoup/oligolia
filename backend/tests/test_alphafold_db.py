"""Tests for backend/services/alphafold_db.py — network calls are always mocked."""

import httpx

from backend.services.alphafold_db import AlphaFoldDBClient

CANNED_ENTRY = [{
    "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P00520-F1-model_v4.pdb",
    "modelCreatedDate": "2022-06-01",
    "globalMetricValue": 91.2,
    "fractionPlddtVeryLow": 0.02,
    "fractionPlddtLow": 0.05,
    "fractionPlddtConfident": 0.20,
    "fractionPlddtVeryHigh": 0.73,
}]
CANNED_PDB_TEXT = "HEADER    TEST\nATOM      1  N   MET A   1      0.0  0.0  0.0  1.00 91.2           N\n"


def test_get_prediction_returns_parsed_result(monkeypatch) -> None:
    def fake_get(self, url, **kwargs):
        if url.endswith("/api/prediction/P00520"):
            return httpx.Response(200, json=CANNED_ENTRY, request=httpx.Request("GET", url))
        if url == CANNED_ENTRY[0]["pdbUrl"]:
            return httpx.Response(200, text=CANNED_PDB_TEXT, request=httpx.Request("GET", url))
        raise AssertionError(f"unexpected URL {url}")

    monkeypatch.setattr(httpx.Client, "get", fake_get)

    prediction = AlphaFoldDBClient().get_prediction("P00520")
    assert prediction is not None
    assert prediction.pdb_text == CANNED_PDB_TEXT
    assert prediction.global_metric_value == 91.2


def test_get_prediction_returns_none_on_404(monkeypatch) -> None:
    def fake_get(self, url, **kwargs):
        return httpx.Response(404, request=httpx.Request("GET", url))

    monkeypatch.setattr(httpx.Client, "get", fake_get)

    assert AlphaFoldDBClient().get_prediction("NOTAREALACCESSION") is None
