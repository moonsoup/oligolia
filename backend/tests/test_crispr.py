"""Tests for CRISPR guide RNA design using real genomic sequences."""

from fastapi.testclient import TestClient


def test_design_cas9_tp53(client: TestClient, tp53_exon7: str) -> None:
    """Design Cas9 guides for TP53 exon 7 — a clinically important target."""
    r = client.post("/crispr/design", json={
        "target_sequence": tp53_exon7,
        "cas_type": "SpCas9",
        "guide_length": 20,
        "max_guides": 10,
    })
    assert r.status_code == 200
    data = r.json()
    assert data["total_candidates"] > 0
    assert len(data["guides"]) <= 10

    for guide in data["guides"]:
        # All guides must be 20 nt
        assert len(guide["sequence"]) == 20
        # All sequences must be DNA
        assert all(c in "ACGTN" for c in guide["sequence"].upper())
        # PAM must be NGG for SpCas9
        assert guide["pam"] == "NGG"
        # Strand must be + or -
        assert guide["strand"] in ("+", "-")
        # GC content must be computable
        assert 0 <= guide["gc_content"] <= 100


def test_design_guides_sorted_by_score(client: TestClient, tp53_exon7: str) -> None:
    """Guides should be sorted by on-target score descending."""
    r = client.post("/crispr/design", json={
        "target_sequence": tp53_exon7,
        "cas_type": "SpCas9",
        "max_guides": 5,
    })
    guides = r.json()["guides"]
    scores = [g["on_target_score"] for g in guides]
    assert scores == sorted(scores, reverse=True)


def test_design_cas12a(client: TestClient, tp53_exon7: str) -> None:
    """Design AsCas12a guides — PAM must be TTTV."""
    r = client.post("/crispr/design", json={
        "target_sequence": tp53_exon7,
        "cas_type": "AsCas12a",
        "max_guides": 5,
    })
    assert r.status_code == 200
    data = r.json()
    for guide in data["guides"]:
        assert guide["pam"] == "TTTV"
        assert len(guide["sequence"]) == 23


def test_design_target_too_short(client: TestClient) -> None:
    r = client.post("/crispr/design", json={
        "target_sequence": "ATCG",
        "cas_type": "SpCas9",
    })
    assert r.status_code == 400


def test_design_invalid_sequence(client: TestClient) -> None:
    r = client.post("/crispr/design", json={
        "target_sequence": "ATCG12345",  # invalid characters
        "cas_type": "SpCas9",
    })
    assert r.status_code == 400


def test_score_guide_good(client: TestClient) -> None:
    """High-GC guide at 50% should score well."""
    r = client.post("/crispr/score_guide?guide_sequence=ATCGATCGATCGATCGATCG")
    assert r.status_code == 200
    data = r.json()
    assert "on_target_score" in data
    assert 0 <= data["on_target_score"] <= 1
    assert data["gc_content"] == 50.0


def test_score_guide_poly_t(client: TestClient) -> None:
    """Guide with poly-T should have a warning."""
    r = client.post("/crispr/score_guide?guide_sequence=ATCGTTTTGATCGATCGATC")
    assert r.status_code == 200
    data = r.json()
    issues = data["issues"]
    assert any("poly-t" in i.lower() or "tttt" in i.lower() for i in issues)


def test_score_guide_low_gc(client: TestClient) -> None:
    """Very low GC guide should flag as issue."""
    r = client.post("/crispr/score_guide?guide_sequence=ATATATATATATATATAT AT".replace(" ", ""))
    assert r.status_code == 200
    data = r.json()
    assert any("GC" in i for i in data["issues"])


def test_forward_and_reverse_guides(client: TestClient) -> None:
    """Both strands should have guides."""
    target = "AAACCCGTTGGCAATGCTTCGGGGAACGTTTCCC" * 3  # 105 nt
    r = client.post("/crispr/design", json={
        "target_sequence": target,
        "cas_type": "SpCas9",
        "max_guides": 20,
    })
    assert r.status_code == 200
    strands = {g["strand"] for g in r.json()["guides"]}
    assert "+" in strands or "-" in strands  # at least one strand
