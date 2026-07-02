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


# ── Gibson assembly (issue #36) ───────────────────────────────────────────────
# Varied 90 bp target so 20 bp overlaps are unambiguous (no accidental repeats).
GIBSON_TARGET = (
    "ATGCGTACGTTGACCTGAGCATGCTAGGCTAACGTTCAGGATCCGATTACGATCGTAG"
    "CTAGCATCGATCGTTAGCACCGGTTAAGCTCAAT"
)


def _split_with_overlaps(seq: str, overlap: int = 20) -> list[str]:
    """Split a circular sequence into 3 fragments sharing `overlap` bp ends."""
    n = len(seq)
    p1, p2 = n // 3, 2 * n // 3
    f0 = seq[0:p1 + overlap]
    f1 = seq[p1:p2 + overlap]
    f2 = seq[p2:n] + seq[0:overlap]
    return [f0, f1, f2]


def test_gibson_assembles_designed_overlaps(client: TestClient) -> None:
    frags = _split_with_overlaps(GIBSON_TARGET, 20)
    # Shuffle input order; assembly must still recover the circle.
    shuffled = [frags[2], frags[0], frags[1]]
    r = client.post("/cloning/gibson", json={"fragments": shuffled, "min_overlap": 15})
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["product"]["is_circular"] is True
    assert len(data["product"]["seq"]) == len(GIBSON_TARGET)
    # Circular product is defined up to rotation.
    assert GIBSON_TARGET in (data["product"]["seq"] + data["product"]["seq"])
    assert sorted(data["order"]) == [0, 1, 2]
    assert len(data["junctions"]) == 3
    assert all(j["overlap_length"] == 20 for j in data["junctions"])


def test_gibson_no_valid_ordering_fails_clearly(client: TestClient) -> None:
    frags = ["A" * 40, "C" * 40, "G" * 40]  # no shared ends
    r = client.post("/cloning/gibson", json={"fragments": frags, "min_overlap": 15})
    assert r.status_code == 400
    assert "No valid circular assembly" in r.json()["detail"]


def test_gibson_ambiguous_ordering_refused(client: TestClient) -> None:
    """Identical linkers on every end -> multiple valid orderings -> refuse."""
    linker = "ACGTACCTGAGTCAGGTTGC"  # 20 bp, non-periodic
    frags = [linker + mid + linker for mid in ("TTTAAACCC", "GAGATCTGA", "TCAGGTACC")]
    r = client.post("/cloning/gibson", json={"fragments": frags, "min_overlap": 15})
    assert r.status_code == 400
    assert "Ambiguous" in r.json()["detail"]


# ── Annotation reindexing onto assembled products (issue #38) ─────────────────

def test_ligation_reindexes_annotations(client: TestClient) -> None:
    """Source annotations shift by each fragment's product offset (hand-checked)."""
    a = {"sequence": "AAAACCCC", "name": "A",
         "left_overhang_type": "blunt", "right_overhang_type": "blunt",
         "annotations": [{"feature_type": "geneA", "start": 0, "end": 4}]}
    b = {"sequence": "GGGGTTTT", "name": "B",
         "left_overhang_type": "blunt", "right_overhang_type": "blunt",
         "annotations": [{"feature_type": "geneB", "start": 2, "end": 6}]}
    r = client.post("/cloning/ligate", json={"fragments": [a, b], "circular": True})
    assert r.status_code == 200, r.text
    anns = {x["feature_type"]: (x["start"], x["end"]) for x in r.json()["product"]["annotations"]}
    assert anns["geneA"] == (0, 4)     # fragment A at offset 0
    assert anns["geneB"] == (10, 14)   # fragment B at offset 8 -> 2+8, 6+8
    assert r.json()["warnings"] == []  # no trimming in ligation


def test_gibson_reindexes_truncates_and_drops(client: TestClient) -> None:
    """Gibson reindex accounts for trimmed overlaps (hand-checked example)."""
    c = GIBSON_TARGET[:40]
    f0seq, f1seq = c[0:30], c[20:40] + c[0:10]  # 2-fragment circle, 10 bp overlaps
    f0 = {"sequence": f0seq, "name": "f0",
          "annotations": [{"feature_type": "gx", "start": 5, "end": 15}]}
    f1 = {"sequence": f1seq, "name": "f1", "annotations": [
        {"feature_type": "drop", "start": 2, "end": 8},    # in trimmed leading overlap
        {"feature_type": "keep", "start": 12, "end": 18},  # retained
        {"feature_type": "span", "start": 8, "end": 14},   # crosses the overlap boundary
    ]}
    r = client.post("/cloning/gibson", json={"fragments": [f0, f1], "min_overlap": 10})
    assert r.status_code == 200, r.text
    # Confirm the designed 10 bp overlaps (guards the hand-checked coordinates).
    assert all(j["overlap_length"] == 10 for j in r.json()["junctions"])
    m = {x["feature_type"]: x for x in r.json()["product"]["annotations"]}
    assert (m["gx"]["start"], m["gx"]["end"]) == (5, 15)     # f0 fully retained
    assert (m["keep"]["start"], m["keep"]["end"]) == (32, 38)
    assert (m["span"]["start"], m["span"]["end"]) == (30, 34)
    assert m["span"]["qualifiers"]["truncated_at_junction"] == "true"
    assert "drop" not in m  # fully inside a trimmed overlap -> dropped
    assert any("drop" in w for w in r.json()["warnings"])


def test_gibson_still_accepts_bare_strings(client: TestClient) -> None:
    """Back-compat: a plain list of sequence strings still assembles."""
    frags = _split_with_overlaps(GIBSON_TARGET, 20)
    r = client.post("/cloning/gibson", json={"fragments": frags, "min_overlap": 15})
    assert r.status_code == 200
    assert r.json()["product"]["annotations"] == []


# ── Golden Gate assembly (issue #37) ──────────────────────────────────────────

def _bsai_part(left_oh: str, body: str, right_oh: str) -> str:
    """A Golden Gate part: two inward BsaI sites releasing left_oh+body+right_oh."""
    return "TTGGTCTCA" + left_oh + body + right_oh + "AGAGACCTT"


def test_golden_gate_assembles_multi_fragment(client: TestClient) -> None:
    p1 = _bsai_part("TCAG", "AAAAAAAA", "AATG")
    p2 = _bsai_part("AATG", "GGGGGGGG", "TCAG")
    r = client.post("/cloning/goldengate", json={"parts": [p1, p2], "enzyme": "BsaI"})
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["product"]["is_circular"] is True
    assert data["product"]["seq"] == "TCAGAAAAAAAA" + "AATGGGGGGGGG"
    assert len(data["junctions"]) == 2
    assert {j["overhang"] for j in data["junctions"]} == {"AATG", "TCAG"}
    assert data["warnings"] == []  # unique, non-palindromic overhangs


def test_golden_gate_palindromic_overhang_flagged(client: TestClient) -> None:
    p1 = _bsai_part("TCAG", "AAAAAAAA", "AATT")   # AATT is palindromic
    p2 = _bsai_part("AATT", "GGGGGGGG", "TCAG")
    r = client.post("/cloning/goldengate", json={"parts": [p1, p2], "enzyme": "BsaI"})
    assert r.status_code == 200, r.text
    assert any("Palindromic" in w and "AATT" in w for w in r.json()["warnings"])


def test_golden_gate_colliding_overhangs_flagged(client: TestClient) -> None:
    p1 = _bsai_part("AATG", "AAAAAAAA", "AATG")   # both junctions AATG -> collision
    p2 = _bsai_part("AATG", "GGGGGGGG", "AATG")
    r = client.post("/cloning/goldengate", json={"parts": [p1, p2], "enzyme": "BsaI"})
    assert r.status_code == 200, r.text
    assert any("Colliding" in w for w in r.json()["warnings"])


def test_golden_gate_incompatible_order_rejected(client: TestClient) -> None:
    p1 = _bsai_part("TCAG", "AAAAAAAA", "AATG")
    p2 = _bsai_part("CCCC", "GGGGGGGG", "TCAG")   # left CCCC != p1 right AATG
    r = client.post("/cloning/goldengate", json={"parts": [p1, p2], "enzyme": "BsaI"})
    assert r.status_code == 400
    assert "Incompatible ends" in r.json()["detail"]


def test_golden_gate_unsupported_enzyme(client: TestClient) -> None:
    p = _bsai_part("TCAG", "AAAAAAAA", "AATG")
    r = client.post("/cloning/goldengate", json={"parts": [p, p], "enzyme": "EcoRI"})
    assert r.status_code == 400
    assert "Type IIS" in r.json()["detail"]


def test_golden_gate_part_without_two_sites_rejected(client: TestClient) -> None:
    bad = "AAAAAAAAAAAAAAAAAAAA"  # no BsaI site -> releases no insert
    good = _bsai_part("TCAG", "AAAAAAAA", "AATG")
    r = client.post("/cloning/goldengate", json={"parts": [bad, good], "enzyme": "BsaI"})
    assert r.status_code == 400
    assert "released" in r.json()["detail"]
