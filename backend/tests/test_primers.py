"""Tests for PCR primer design and restriction enzyme analysis using real sequences."""

from fastapi.testclient import TestClient


def test_design_primers_hbb(client: TestClient, test_template: str) -> None:
    """Design primers for the test template (contains common restriction sites)."""
    r = client.post("/primers/design", json={
        "template": test_template,
        "product_min": 80,
        "product_max": 250,
        "primer_len_min": 18,
        "primer_len_max": 22,
        "tm_min": 50.0,
        "tm_max": 70.0,
        "max_pairs": 5,
    })
    assert r.status_code == 200
    pairs = r.json()
    assert len(pairs) >= 1
    for pair in pairs:
        assert pair["forward"]["direction"] == "forward"
        assert pair["reverse"]["direction"] == "reverse"
        assert pair["product_size"] >= 80
        assert pair["product_size"] <= 250
        fwd_len = len(pair["forward"]["sequence"])
        assert 18 <= fwd_len <= 22
        # Tm should be within requested range
        assert 50.0 <= pair["forward"]["tm"] <= 70.0
        assert 50.0 <= pair["reverse"]["tm"] <= 70.0


def test_primer_gc_content_range(client: TestClient, test_template: str) -> None:
    r = client.post("/primers/design", json={
        "template": test_template,
        "product_min": 80,
        "product_max": 300,
        "gc_min": 40.0,
        "gc_max": 65.0,
        "max_pairs": 3,
    })
    assert r.status_code == 200
    for pair in r.json():
        assert 40.0 <= pair["forward"]["gc_content"] <= 65.0
        assert 40.0 <= pair["reverse"]["gc_content"] <= 65.0


def test_primers_sorted_by_penalty(client: TestClient, test_template: str) -> None:
    r = client.post("/primers/design", json={
        "template": test_template,
        "product_min": 80,
        "product_max": 300,
        "max_pairs": 5,
    })
    pairs = r.json()
    penalties = [p["penalty"] for p in pairs]
    assert penalties == sorted(penalties)


def test_template_too_short(client: TestClient) -> None:
    r = client.post("/primers/design", json={
        "template": "ATCGATCG",
        "product_min": 100,
    })
    assert r.status_code == 400


def test_restriction_sites_ecori(client: TestClient, test_template: str) -> None:
    """Template starts with AATTC — EcoRI site (GAATTC) should be found."""
    r = client.post("/primers/restriction_sites", json={"template": test_template, "enzymes": ["EcoRI"]})
    assert r.status_code == 200
    sites = r.json()
    ecori = next((s for s in sites if s["enzyme"] == "EcoRI"), None)
    assert ecori is not None
    assert ecori["count"] >= 1


def test_restriction_sites_bamhi(client: TestClient, test_template: str) -> None:
    """Template contains GGATCC (BamHI) multiple times."""
    r = client.post("/primers/restriction_sites", json={"template": test_template})
    assert r.status_code == 200
    sites = r.json()
    bamhi = next((s for s in sites if s["enzyme"] == "BamHI"), None)
    assert bamhi is not None
    assert bamhi["count"] >= 2  # template has GGATCC twice


def test_restriction_sites_all_enzymes(client: TestClient, test_template: str) -> None:
    """Should return only enzymes that actually cut the template."""
    r = client.post("/primers/restriction_sites", json={"template": test_template})
    assert r.status_code == 200
    sites = r.json()
    assert len(sites) >= 1
    for site in sites:
        assert site["count"] >= 1
        assert len(site["positions"]) == site["count"]


def test_restriction_sites_empty_template(client: TestClient) -> None:
    """Edge case: short template with no enzyme sites."""
    r = client.post("/primers/restriction_sites", json={"template": "AAAAAAAAAA"})
    assert r.status_code == 200
    assert r.json() == []
