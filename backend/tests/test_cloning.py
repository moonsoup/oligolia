"""Tests for restriction-ligation assembly (issue #35).

Fragments are produced by the real digest endpoint (issue #34) and fed into
the ligation endpoint, exercising the two features end-to-end.
"""

from fastapi.testclient import TestClient

from ..routers.primers import DigestRequest, digest

# Circular plasmid with exactly two EcoRI sites -> circular digest gives 2 fragments.
TWO_ECORI = "GAATTC" + "A" * 12 + "GAATTC" + "T" * 12  # 36 bp


def _frag_dicts_in_order(template: str, enzymes: list[str], circular: bool) -> list[dict]:
    res = digest(DigestRequest(template=template, enzymes=enzymes, is_circular=circular))
    frags = sorted(res.fragments, key=lambda f: f.start)  # circular order
    return [f.model_dump() for f in frags]


def test_religate_circular_digest_restores_the_circle(client: TestClient) -> None:
    """Cutting a circular plasmid then ligating the pieces rebuilds it."""
    frags = _frag_dicts_in_order(TWO_ECORI, ["EcoRI"], circular=True)
    assert len(frags) == 2
    r = client.post("/cloning/ligate", json={"fragments": frags, "circular": True})
    assert r.status_code == 200, r.text
    product = r.json()["product"]
    assert product["is_circular"] is True
    assert len(product["seq"]) == len(TWO_ECORI)
    # A circular product is defined up to rotation — the original must appear
    # in the doubled product sequence.
    assert TWO_ECORI in (product["seq"] + product["seq"])
    # Two sticky 5' AATT junctions closed the circle.
    assert len(r.json()["junctions"]) == 2
    assert all(j["kind"] == "sticky-5'" and j["overhang"] == "AATT"
               for j in r.json()["junctions"])


def test_incompatible_ends_are_rejected(client: TestClient) -> None:
    """An EcoRI end cannot ligate to a BamHI end — clear error, no silent join."""
    ecori = {"sequence": "AATTCGGG", "name": "ecoRI-frag",
             "left_overhang": "AATT", "left_overhang_type": "5'",
             "right_overhang": "AATT", "right_overhang_type": "5'"}
    bamhi = {"sequence": "GATCCTTT", "name": "bamHI-frag",
             "left_overhang": "GATC", "left_overhang_type": "5'",
             "right_overhang": "GATC", "right_overhang_type": "5'"}
    r = client.post("/cloning/ligate", json={"fragments": [ecori, bamhi], "circular": True})
    assert r.status_code == 400
    assert "Incompatible ends" in r.json()["detail"]
    assert "AATT" in r.json()["detail"] and "GATC" in r.json()["detail"]


def test_blunt_fragments_ligate(client: TestClient) -> None:
    """Blunt ends (e.g. SmaI) ligate to each other."""
    a = {"sequence": "AAACCC", "name": "a",
         "left_overhang_type": "blunt", "right_overhang_type": "blunt"}
    b = {"sequence": "GGGTTT", "name": "b",
         "left_overhang_type": "blunt", "right_overhang_type": "blunt"}
    r = client.post("/cloning/ligate", json={"fragments": [a, b], "circular": True})
    assert r.status_code == 200
    assert r.json()["product"]["seq"] == "AAACCCGGGTTT"
    assert all(j["kind"] == "blunt" for j in r.json()["junctions"])


def test_linear_ligation_skips_closing_junction(client: TestClient) -> None:
    """circular=False checks only the internal junction, not last->first."""
    frags = _frag_dicts_in_order(TWO_ECORI, ["EcoRI"], circular=True)
    r = client.post("/cloning/ligate", json={"fragments": frags, "circular": False})
    assert r.status_code == 200
    assert r.json()["product"]["is_circular"] is False
    assert len(r.json()["junctions"]) == 1  # one internal junction, no closure


def test_at_least_two_fragments_required(client: TestClient) -> None:
    frag = {"sequence": "AATTCGGG",
            "left_overhang": "AATT", "left_overhang_type": "5'",
            "right_overhang": "AATT", "right_overhang_type": "5'"}
    r = client.post("/cloning/ligate", json={"fragments": [frag], "circular": True})
    assert r.status_code == 422  # pydantic min_length


def test_free_linear_ends_cannot_ligate(client: TestClient) -> None:
    """A 'none' (free terminus) end forms no defined junction."""
    a = {"sequence": "AAAA", "name": "a",
         "left_overhang_type": "none", "right_overhang_type": "none"}
    b = {"sequence": "TTTT", "name": "b",
         "left_overhang_type": "none", "right_overhang_type": "none"}
    r = client.post("/cloning/ligate", json={"fragments": [a, b], "circular": True})
    assert r.status_code == 400
