"""Tests for backend/services/esmfold.py — network calls are always mocked."""

import httpx
import pytest

from backend.services.esmfold import ESMFoldClient, MAX_LENGTH


def test_predict_returns_pdb_text(monkeypatch) -> None:
    canned_pdb = "HEADER    TEST\nATOM      1  N   MET A   1      0.0  0.0  0.0  1.00 50.0           N\n"

    def fake_post(self, url, content=None, **kwargs):
        assert content == "MKV"
        return httpx.Response(200, text=canned_pdb, request=httpx.Request("POST", url))

    monkeypatch.setattr(httpx.Client, "post", fake_post)

    result = ESMFoldClient().predict("mkv")
    assert result == canned_pdb


def test_predict_rejects_over_length_sequence_without_network_call(monkeypatch) -> None:
    def fail_post(self, *args, **kwargs):
        raise AssertionError("should not make a network call for an over-length sequence")

    monkeypatch.setattr(httpx.Client, "post", fail_post)

    too_long = "A" * (MAX_LENGTH + 1)
    with pytest.raises(ValueError, match=r"400"):
        ESMFoldClient().predict(too_long)
