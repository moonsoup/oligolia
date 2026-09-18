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


def test_off_targets_disabled_by_default(client: TestClient, tp53_exon7: str) -> None:
    """Without check_off_targets, guides carry no off-target fields."""
    r = client.post("/crispr/design", json={
        "target_sequence": tp53_exon7,
        "cas_type": "SpCas9",
        "max_guides": 5,
    })
    assert r.status_code == 200
    for g in r.json()["guides"]:
        assert g["off_target_count"] is None
        assert g["off_target_summary"] is None
        assert g["specificity_score"] is None


def test_off_targets_populated(client: TestClient, tp53_exon7: str) -> None:
    """check_off_targets fills mismatch buckets and a specificity score."""
    r = client.post("/crispr/design", json={
        "target_sequence": tp53_exon7,
        "cas_type": "SpCas9",
        "max_guides": 5,
        "check_off_targets": True,
    })
    assert r.status_code == 200
    guides = r.json()["guides"]
    assert guides
    for g in guides:
        assert g["off_target_count"] is not None
        summary = g["off_target_summary"]
        assert set(summary.keys()) == {"0", "1", "2", "3"}
        assert g["off_target_count"] == sum(summary.values())
        assert 0 <= g["specificity_score"] <= 100


def test_off_targets_find_paralog(client: TestClient) -> None:
    """An exact duplicate site in a reference sequence is flagged as off-target."""
    guide = "ACGTACGTACGTACGTACGT"
    # A site is protospacer + N + GG: the N belongs to the PAM, not the guide (#53).
    target = "GGGG" + guide + "AGG" + "CCCC"           # on-target site (protospacer + NGG)
    paralog = "AAAA" + guide + "AGG" + "TTTT"          # identical protospacer elsewhere
    r = client.post("/crispr/design", json={
        "target_sequence": target,
        "cas_type": "SpCas9",
        "max_guides": 20,
        "check_off_targets": True,
        "reference_sequences": [paralog],
    })
    assert r.status_code == 200
    hit = next((g for g in r.json()["guides"] if g["sequence"] == guide), None)
    assert hit is not None
    assert hit["off_target_summary"]["0"] >= 1  # the paralog's exact match
    assert hit["specificity_score"] < 100


def test_off_targets_skipped_for_cas13(client: TestClient, tp53_exon7: str) -> None:
    """Cas13 targets RNA; genomic off-target scan does not apply."""
    r = client.post("/crispr/design", json={
        "target_sequence": tp53_exon7,
        "cas_type": "LwaCas13a",
        "max_guides": 5,
        "check_off_targets": True,
    })
    assert r.status_code == 200
    for g in r.json()["guides"]:
        assert g["off_target_count"] is None
        assert g["specificity_score"] is None


def test_spcas9_guide_is_the_protospacer_and_excludes_the_pam(client: TestClient) -> None:
    """#53. The SpCas9 PAM is N-GG, and `(?=(.{20})GG)` captured the N as the
    guide's last base — so every returned guide was shifted one base along the
    target and carried a PAM base, and the real protospacer was never returned
    at all. The mismatch lands in the PAM-proximal seed, where SpCas9 is least
    tolerant, so ordered guides largely would not cut.
    """
    protospacer = "GACGTTACGATCGGATCCAT"
    target = "CCCCC" + protospacer + "TGG" + "CCCCC"

    r = client.post("/crispr/design", json={
        "target_sequence": target, "cas_type": "SpCas9", "max_guides": 20})
    assert r.status_code == 200
    guides = r.json()["guides"]
    found = [g for g in guides if g["strand"] == "+"]

    assert protospacer in [g["sequence"] for g in found], [g["sequence"] for g in found]
    # the shifted read: the last base of the guide taken from the PAM
    assert "ACGTTACGATCGGATCCATT" not in [g["sequence"] for g in guides]
    hit = next(g for g in found if g["sequence"] == protospacer)
    assert hit["position"] == 5, hit


def test_spcas9_reverse_strand_guide_is_also_the_protospacer(client: TestClient) -> None:
    """The reverse-strand branch carried the same off-by-one (#53)."""
    protospacer = "GACGTTACGATCGGATCCAT"
    # on the minus strand: the reverse complement of (protospacer + N + GG)
    site = protospacer + "TGG"
    rc_site = site.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    target = "AAAAA" + rc_site + "AAAAA"

    r = client.post("/crispr/design", json={
        "target_sequence": target, "cas_type": "SpCas9", "max_guides": 20})
    assert r.status_code == 200
    minus = [g["sequence"] for g in r.json()["guides"] if g["strand"] == "-"]
    assert protospacer in minus, minus


# --- #64: non-Cas9 nuclease handling ---

def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_cas13_includes_the_final_window(client: TestClient) -> None:
    """#64.1: range(0, len - 22) stopped one short, so a 24 nt target gave 2 of 3."""
    target = "ACGT" * 6  # 24 nt
    r = client.post("/crispr/design", json={
        "target_sequence": target, "cas_type": "LwaCas13a", "max_guides": 50,
    })
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["total_candidates"] == 3, data["total_candidates"]
    positions = sorted(g["position"] for g in data["guides"])
    assert positions == [0, 1, 2]


def test_cas13_window_count_is_general(client: TestClient) -> None:
    for n in (22, 23, 30, 45):
        target = ("ACGTGGCATC" * 10)[:n]
        r = client.post("/crispr/design", json={
            "target_sequence": target, "cas_type": "LwaCas13a", "max_guides": 50,
        })
        assert r.status_code == 200, r.text
        assert r.json()["total_candidates"] == n - 22 + 1, (n, r.json()["total_candidates"])


def test_cas12a_finds_minus_strand_sites(client: TestClient) -> None:
    """#64.2: only the + strand was scanned, so a minus-strand-only target gave 0."""
    # Build a target whose TTTV PAM + 23 nt guide exists only on the minus strand.
    minus = "TTTA" + "ACGTGGCATCACGTGGCATCACG"  # PAM + 23 nt, on the reverse strand
    target = _rc(minus)
    assert "TTT" not in target, "the forward strand must carry no TTTV site for this test"

    r = client.post("/crispr/design", json={
        "target_sequence": target, "cas_type": "AsCas12a", "max_guides": 50,
    })
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["total_candidates"] >= 1, "no minus-strand Cas12a site found (#64)"
    assert any(g["strand"] == "-" for g in data["guides"])


def test_guide_length_is_honoured_when_given(client: TestClient) -> None:
    """#64.3: the parameter had no effect for any Cas type."""
    target = "ACGTGGCATC" * 10
    for length in (17, 18, 20, 24):
        r = client.post("/crispr/design", json={
            "target_sequence": target, "cas_type": "SpCas9",
            "guide_length": length, "max_guides": 50,
        })
        assert r.status_code == 200, r.text
        for g in r.json()["guides"]:
            assert len(g["sequence"]) == length, (length, g["sequence"])


def test_guide_length_defaults_stay_canonical_per_nuclease(client: TestClient) -> None:
    """Not passing guide_length must keep each nuclease's usual length."""
    target = "TTTA" + "ACGTGGCATC" * 10
    r = client.post("/crispr/design", json={
        "target_sequence": target, "cas_type": "AsCas12a", "max_guides": 50,
    })
    assert r.status_code == 200, r.text
    for g in r.json()["guides"]:
        assert len(g["sequence"]) == 23, g["sequence"]

    r = client.post("/crispr/design", json={
        "target_sequence": target, "cas_type": "LwaCas13a", "max_guides": 50,
    })
    assert r.status_code == 200, r.text
    for g in r.json()["guides"]:
        assert len(g["sequence"]) == 22, g["sequence"]


def test_every_guide_re_derives_from_the_raw_target(client: TestClient) -> None:
    """The invariant that catches this whole defect class, including #53.

    For every guide of every Cas type: re-extract it from the raw target by its
    reported position and strand, and check the PAM grammar around it. A shifted
    guide, a guide that includes a PAM base, or a reverse-strand index that maps
    the wrong way all fail here.
    """
    target = (
        "TTTAACGTGGCATCACGTGGCATCACGATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGG"
        "ACGTGGCATCACGTGGCATCACGTTTCGGCCATGGACGTAGGCATCACGTGGCATCACGTGG"
    )
    for cas, pam_check in (
        ("SpCas9", "NGG"),
        ("SpCas9-HF1", "NGG"),
        ("AsCas12a", "TTTV"),
        ("LwaCas13a", None),
    ):
        r = client.post("/crispr/design", json={
            "target_sequence": target, "cas_type": cas, "max_guides": 50,
        })
        assert r.status_code == 200, (cas, r.text)
        guides = r.json()["guides"]
        assert guides, f"{cas} returned no guides on this target"

        for g in guides:
            pos, strand, seq = g["position"], g["strand"], g["sequence"]
            L = len(seq)
            assert 0 <= pos and pos + L <= len(target), (cas, g)

            window = target[pos:pos + L]
            expected = window if strand == "+" else _rc(window)
            assert expected == seq, f"{cas} {strand} guide does not re-derive: {g}"

            if pam_check == "NGG":
                # PAM sits 3' of the protospacer on that strand.
                if strand == "+":
                    pam = target[pos + L:pos + L + 3]
                else:
                    pam = _rc(target[pos - 3:pos])
                assert len(pam) == 3, (cas, g, pam)
                assert pam[1:] == "GG", f"{cas} {strand} PAM is not NGG: {pam!r} for {g}"
                assert not seq.endswith(pam), "the guide must not include the PAM (#53)"
            elif pam_check == "TTTV":
                if strand == "+":
                    pam = target[pos - 4:pos]
                else:
                    pam = _rc(target[pos + L:pos + L + 4])
                assert len(pam) == 4, (cas, g, pam)
                assert pam[:3] == "TTT" and pam[3] in "ACG", f"{cas} PAM is not TTTV: {pam!r}"
