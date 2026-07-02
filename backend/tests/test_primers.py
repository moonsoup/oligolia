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


# ── Circular topology (issue #21) ─────────────────────────────────────────────
# Synthetic 32 bp construct: an internal EcoRI site (GAATTC at index 11) plus an
# origin-spanning one — the template ends with "GAATT" and starts with "C", so
# reading around the junction spells GAATTC at real start index 27.
CIRCULAR_ECORI = "C" + "A" * 10 + "GAATTC" + "A" * 10 + "GAATT"  # len 32


def test_circular_finds_origin_spanning_site(client: TestClient) -> None:
    linear = client.post("/primers/restriction_sites",
                         json={"template": CIRCULAR_ECORI, "enzymes": ["EcoRI"]}).json()
    circular = client.post("/primers/restriction_sites",
                           json={"template": CIRCULAR_ECORI, "enzymes": ["EcoRI"],
                                 "is_circular": True}).json()
    lin_pos = next(s for s in linear if s["enzyme"] == "EcoRI")["positions"]
    circ_pos = next(s for s in circular if s["enzyme"] == "EcoRI")["positions"]
    assert lin_pos == [11]              # linear scan misses the wrap
    assert circ_pos == [11, 27]         # circular scan finds the origin-spanning site


def test_circular_digest_n_fragments_for_n_cuts(client: TestClient) -> None:
    """A circular molecule with N cuts yields N fragments that tile the circle."""
    r = client.post("/primers/digest",
                    json={"template": CIRCULAR_ECORI, "enzymes": ["EcoRI"], "is_circular": True})
    assert r.status_code == 200
    data = r.json()
    # Real EcoRI geometry (G^AATTC): internal cut at 12 and origin-spanning at 28.
    assert sorted(data["cut_positions"]) == [12, 28]
    frags = data["fragments"]
    assert len(frags) == 2  # N cuts -> N fragments
    assert sum(f["length"] for f in frags) == len(CIRCULAR_ECORI)
    # The origin-spanning fragment reconstructs across the junction.
    wrap = next(f for f in frags if f["end"] < f["start"])
    assert wrap["sequence"] == CIRCULAR_ECORI[wrap["start"]:] + CIRCULAR_ECORI[:wrap["end"]]


def test_linear_digest_unchanged(client: TestClient) -> None:
    """Default (linear) behavior: N cuts -> N+1 fragments."""
    r = client.post("/primers/digest",
                    json={"template": CIRCULAR_ECORI, "enzymes": ["EcoRI"]})
    data = r.json()
    assert data["cut_positions"] == [12]        # only the internal site cuts, real geometry
    assert len(data["fragments"]) == 2          # N+1 = 2 fragments for 1 cut
    assert sum(f["length"] for f in data["fragments"]) == len(CIRCULAR_ECORI)


# ── Real cut geometry + overhangs (issue #34) ─────────────────────────────────

def test_digest_real_cut_geometry_ecori(client: TestClient) -> None:
    """EcoRI cuts G^AATTC, not after the whole site — verified vs NEB geometry."""
    r = client.post("/primers/digest",
                    json={"template": "AAAGAATTCTTT", "enzymes": ["EcoRI"]})
    data = r.json()
    assert data["cut_positions"] == [4]  # after the G in AAAG^AATTC
    seqs = {f["sequence"] for f in data["fragments"]}
    assert seqs == {"AAAG", "AATTCTTT"}


def test_digest_5prime_overhang(client: TestClient) -> None:
    r = client.post("/primers/digest",
                    json={"template": "AAAGAATTCTTT", "enzymes": ["EcoRI"]})
    frags = r.json()["fragments"]
    # Both ends flanking the cut carry the 5' AATT overhang.
    left = next(f for f in frags if f["start"] == 0)
    right = next(f for f in frags if f["start"] == 4)
    assert left["right_overhang_type"] == "5'" and left["right_overhang"] == "AATT"
    assert right["left_overhang_type"] == "5'" and right["left_overhang"] == "AATT"
    # Free linear termini have no overhang.
    assert left["left_overhang_type"] == "none"
    assert right["right_overhang_type"] == "none"


def test_digest_3prime_overhang(client: TestClient) -> None:
    r = client.post("/primers/digest",
                    json={"template": "AAACTGCAGTTT", "enzymes": ["PstI"]})
    frags = r.json()["fragments"]
    up = next(f for f in frags if f["start"] == 0)
    assert up["right_overhang_type"] == "3'" and up["right_overhang"] == "TGCA"


def test_digest_blunt_cutter(client: TestClient) -> None:
    r = client.post("/primers/digest",
                    json={"template": "AAACCCGGGTTT", "enzymes": ["SmaI"]})
    frags = r.json()["fragments"]
    inner = [f for f in frags if f["left_overhang_type"] == "blunt" or f["right_overhang_type"] == "blunt"]
    assert inner  # the cut produced blunt ends
    for f in frags:
        assert f["left_overhang"] == "" and f["right_overhang"] == ""
