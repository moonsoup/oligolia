"""#62: primer filter defects, separate from #50 (Tm), #51 and #52.

  1. the hairpin filter flagged 56% of random 20-mers, including
     GAATTCAAAAAAAAAAAAAA — GAATTC is self-complementary, but its two halves are
     ADJACENT, so there is no loop and it cannot fold back on itself. The check
     had no loop requirement at all;
  2. reverse candidates skipped the 3'-GC-clamp and homopolymer checks that
     forward candidates got, so 14 of 20 returned pairs had reverse primers
     violating the forward rules;
  3. pairing was O(F x R) — 5.2 s at 3 kb.
"""

from __future__ import annotations

import random
import time

from fastapi.testclient import TestClient

from backend.routers.primers import _acceptable_primer, _has_hairpin


def _random_20mers(n: int, seed: int = 7) -> list[str]:
    rng = random.Random(seed)
    return ["".join(rng.choice("ACGT") for _ in range(20)) for _ in range(n)]


# ── the hairpin filter ───────────────────────────────────────────────────────

def test_an_adjacent_palindrome_is_not_a_hairpin() -> None:
    """The issue's example. GAATTC is self-complementary and cannot fold: its
    two arms touch, so there is no loop."""
    assert _has_hairpin("GAATTCAAAAAAAAAAAAAA") is False


def test_a_real_stem_and_loop_is_a_hairpin() -> None:
    """A 5-bp arm, a loop of 5, then the arm's reverse complement."""
    assert _has_hairpin("GCGCA" + "AAAAA" + "TGCGC") is True


def test_a_four_base_stem_is_not_enough() -> None:
    """Deliberate: MIN_HAIRPIN_STEM is 5.

    A 4-mer has only 256 possibilities, so ordinary sequence pairs up somewhere
    by chance. Measured on 2000 random 20-mers: a 4-bp stem flags 16.6%, a 5-bp
    stem 2.9%, a 6-bp stem 0.5%. Rejecting one primer in six for a 4-bp stem is
    not worth it — real tools score hairpins by free energy, and this is a length
    heuristic set where it does not reject most of the template.
    """
    assert _has_hairpin("GCGC" + "AAAAA" + "GCGC") is False


def test_a_stem_with_too_short_a_loop_is_not_a_hairpin() -> None:
    """A 1-nt loop is not physically foldable."""
    assert _has_hairpin("GCGCA" + "A" + "TGCGC") is False


def test_a_clean_primer_is_not_a_hairpin() -> None:
    assert _has_hairpin("ACGTGGCATCACGATGGCCT") is False


def test_the_rejection_rate_on_random_20mers_is_not_half() -> None:
    """The measurable claim in the issue: 56% of random 20-mers were flagged."""
    seqs = _random_20mers(500)
    flagged = sum(1 for s in seqs if _has_hairpin(s))
    rate = flagged / len(seqs)
    assert rate < 0.20, f"{rate:.0%} of random 20-mers flagged as hairpins"


def test_a_longer_stem_still_counts() -> None:
    assert _has_hairpin("GGCCTT" + "AAGTCA" + "AAGGCC") is True


# ── one filter for both orientations ────────────────────────────────────────

def test_the_shared_filter_rejects_a_3prime_gc_clamp_over_two() -> None:
    assert _acceptable_primer("ACGTGGCATCACGATGGCCG") is False   # 3 of last 3 are G/C
    assert _acceptable_primer("ACGTGGCATCACGATGGCCT") is True    # 2 of last 3


def test_the_shared_filter_rejects_homopolymer_runs() -> None:
    for bad in ("ACGTGGCATAAAAATGGCCT", "ACGTGGCATTTTTTTGGCCT",
                "ACGTGGCATGGGGATGGCCT", "ACGTGGCATCCCCATGGCCT"):
        assert _acceptable_primer(bad) is False, bad


def test_the_shared_filter_rejects_a_hairpin() -> None:
    assert _acceptable_primer("GCGCA" + "AGTCA" + "TGCGCTAGT") is False


def test_the_reverse_primer_of_the_issue_is_rejected() -> None:
    """GCTGCCCCTACGGATCGCA — CCCC run, and a 3-of-3 GC clamp."""
    assert _acceptable_primer("GCTGCCCCTACGGATCGCA") is False


# ── end to end ───────────────────────────────────────────────────────────────

TEMPLATE = (
    "ATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACGTGGCATCACGTTT"
    "CGGCCATGGACGTAGGCATCACGTGGCATCACGTGGCCGGCGGCGGCGCTGCTGCAAGCTT"
    "ACGTGGCATCACGATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACG"
    "TTCGGCCATGGACGTAGGCATCACGTGGCATCACGTTTAGGCCATGGACGTGGCATCACGA"
)


def _pairs(client: TestClient, template: str, **params) -> list[dict]:
    r = client.post("/primers/design", json={
        "template": template, "tm_min": 45.0, "tm_max": 75.0,
        "product_min": 80, "product_max": 200, "max_pairs": 20, **params,
    })
    assert r.status_code == 200, r.text
    return r.json()


def test_every_returned_reverse_primer_obeys_the_forward_rules(client: TestClient) -> None:
    """#62.2: 14 of 20 pairs had reverse primers that violated them."""
    pairs = _pairs(client, TEMPLATE)
    assert pairs, "no pairs returned for this template"

    offenders = [
        p["reverse"]["sequence"] for p in pairs
        if not _acceptable_primer(p["reverse"]["sequence"])
    ]
    assert not offenders, f"{len(offenders)} of {len(pairs)} reverse primers break the rules: {offenders[:5]}"


def test_every_returned_forward_primer_obeys_them_too(client: TestClient) -> None:
    pairs = _pairs(client, TEMPLATE)
    offenders = [
        p["forward"]["sequence"] for p in pairs
        if not _acceptable_primer(p["forward"]["sequence"])
    ]
    assert not offenders, offenders[:5]


def test_a_3kb_template_completes_promptly(client: TestClient) -> None:
    """#62.3: 5.2 s at 3 kb, and worse as product_max grows."""
    template = (TEMPLATE * 13)[:3000]
    started = time.monotonic()
    _pairs(client, template, product_max=1000)
    elapsed = time.monotonic() - started
    assert elapsed < 3.0, f"3 kb design took {elapsed:.1f}s"


def test_pairs_still_respect_the_product_window(client: TestClient) -> None:
    """Guard the optimisation: it must not change which pairs are legal."""
    for lo, hi in ((80, 120), (150, 200), (100, 300)):
        for p in _pairs(client, TEMPLATE, product_min=lo, product_max=hi):
            assert lo <= p["product_size"] <= hi, (lo, hi, p["product_size"])
            assert p["reverse"]["position"] > p["forward"]["position"]


def test_the_reverse_primer_is_the_reverse_complement_of_the_template(client: TestClient) -> None:
    """Re-derivation invariant, as for CRISPR guides and repeats."""
    def rc(s: str) -> str:
        return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    for p in _pairs(client, TEMPLATE):
        rev = p["reverse"]
        window = TEMPLATE[rev["position"]:rev["position"] + rev["length"]]
        assert rev["sequence"] == rc(window), rev
        fwd = p["forward"]
        assert fwd["sequence"] == TEMPLATE[fwd["position"]:fwd["position"] + fwd["length"]], fwd


# --- the pairing optimisation must not change which pairs come back ---

def test_the_pruned_search_returns_exactly_the_brute_force_top_pairs(client: TestClient) -> None:
    """#62.3's fix is an EXACT prune, and this proves it rather than asserting it.

    `penalty = tm_diff + 0.1 * gc_diff` and gc_diff >= 0, so penalty >= tm_diff
    always. Once the result heap is full, a pair whose tm_diff already exceeds the
    worst kept penalty cannot beat it, so it can be skipped without computing the
    GC term. Nothing heuristic about it — but an "optimisation" that quietly
    changes results is the defect class this project keeps finding, so it is
    checked against a full enumeration.

    Measured on this template: 1,172,601 legal pairs, identical penalties and
    identical pair set. 3 kb design went from 83 s to 2.9 s.
    """
    from backend.routers import primers as P

    template = (TEMPLATE * 3)[:750]
    req = P.PrimerDesignRequest(
        template=template, tm_min=45.0, tm_max=75.0,
        product_min=80, product_max=300, max_pairs=20,
    )
    got = P.design_primers(req)
    assert got, "no pairs returned"

    # Independent enumeration, same filters, no pruning at all.
    tpl = template.upper()
    fwd, rev = [], []
    for length in range(req.primer_len_min, req.primer_len_max + 1):
        for pos in range(0, len(tpl) - length + 1):
            seq = tpl[pos:pos + length]
            gc = P._gc(seq)
            if not (req.gc_min <= gc <= req.gc_max) or not P._acceptable_primer(seq):
                continue
            tm = P._tm_nearest_neighbor(seq)
            if req.tm_min <= tm <= req.tm_max:
                fwd.append((seq, pos, length, round(tm, 1), round(gc, 1)))
        for pos in range(length, len(tpl) + 1):
            seq = P._reverse_complement(tpl[pos - length:pos])
            gc = P._gc(seq)
            if not (req.gc_min <= gc <= req.gc_max) or not P._acceptable_primer(seq):
                continue
            tm = P._tm_nearest_neighbor(seq)
            if req.tm_min <= tm <= req.tm_max:
                rev.append((seq, pos - length, length, round(tm, 1), round(gc, 1)))

    every = []
    for fs, fp, fl, ftm, fgc in fwd:
        for rs, rp, rl, rtm, rgc in rev:
            product = rp + rl - fp
            if not (req.product_min <= product <= req.product_max) or rp <= fp:
                continue
            tm_diff = abs(ftm - rtm)
            if tm_diff > 5:
                continue
            every.append((round(tm_diff + abs(fgc - rgc) * 0.1, 3), fs, rs, product))
    every.sort(key=lambda x: x[0])
    assert len(every) > 1000, f"only {len(every)} legal pairs — not a real test"

    mine = [(round(g.penalty, 3), g.forward.sequence, g.reverse.sequence, g.product_size)
            for g in got]
    assert [p[0] for p in every[:len(got)]] == [m[0] for m in mine], "penalties differ"
    assert sorted(every[:len(got)]) == sorted(mine), "the pruned search returned different pairs"
