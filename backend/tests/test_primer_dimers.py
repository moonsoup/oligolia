"""#52: no forward/reverse 3'-dimer check.

A 3' primer-dimer is a property of the PAIR, not of either primer, so no amount
of per-primer filtering catches it — and it is one of the commonest reasons a PCR
consumes itself without amplifying the target. `design_primers` checked each
primer for a self-hairpin and never compared the two to each other.

The check is pure, so these run offline.
"""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from backend.routers.primers import (
    MAX_3PRIME_DIMER,
    forms_3prime_dimer,
    three_prime_dimer_length,
)


def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# ── the pure check ───────────────────────────────────────────────────────────

def test_two_primers_with_complementary_3_prime_ends_form_a_dimer() -> None:
    """The textbook case: the last bases of each anneal to the other."""
    fwd = "ACGTGGCATCACGATGGCCT"
    # Make the reverse primer's 3' tail the reverse complement of the forward's.
    rev = "TTACCATGCATTTA" + _rc(fwd[-6:])
    assert forms_3prime_dimer(fwd, rev) is True
    assert three_prime_dimer_length(fwd, rev) >= 6


def test_a_non_complementary_pair_does_not() -> None:
    """A single complementary base is length 1 and is correctly below threshold.

    Any two primers whose last bases happen to be A/T or G/C give an overlap of
    1 — that is real and harmless, which is why the threshold exists rather than
    a zero test.
    """
    fwd, rev = "ACGTGGCATCACGATGGCCT", "TTTATTTATTTATTTATTTA"
    assert forms_3prime_dimer(fwd, rev) is False
    assert three_prime_dimer_length(fwd, rev) <= MAX_3PRIME_DIMER


def test_the_overlap_length_is_reported() -> None:
    fwd = "ACGTGGCATCACGATGGCCTAG"
    for n in (3, 4, 5, 6, 7):
        rev = "TTATTATTATTATT" + _rc(fwd[-n:])
        assert three_prime_dimer_length(fwd, rev) >= n, (n, rev)


def test_only_the_3_prime_ends_count() -> None:
    """Complementarity at the 5' ends does not prime; the 3' end is what extends.

    Both primers here begin with a fully complementary 10-mer and end with the
    same non-complementary tail, so the overlap is 0 despite the 5' match.
    """
    fwd = "GGCCTTACGT" + "CCCATAGTCAG"
    rev = _rc("GGCCTTACGT") + "GGGATAGTCAG"
    assert three_prime_dimer_length(fwd, rev) == 0
    assert forms_3prime_dimer(fwd, rev) is False


def test_the_threshold_is_respected() -> None:
    """Short chance overlaps must not disqualify every pair."""
    fwd = "ACGTGGCATCACGATGGCCT"
    for n in (1, 2, 3):
        short = "TTATTATTATTATTATT" + _rc(fwd[-n:])
        assert forms_3prime_dimer(fwd, short) is False, (n, short)
    # One more base and it is a finding.
    assert forms_3prime_dimer(fwd, "TTATTATTATTATTATT" + _rc(fwd[-4:])) is True


def test_a_primer_against_itself_is_a_self_dimer() -> None:
    """Same check, both arguments the same primer — the self-dimer case."""
    selfy = "ACGTACGTAC" + _rc("ACGTACGTAC")
    assert forms_3prime_dimer(selfy, selfy) is True


# ── end to end ───────────────────────────────────────────────────────────────

TEMPLATE = (
    "ATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACGTGGCATCACGTTT"
    "CGGCCATGGACGTAGGCATCACGTGGCATCACGTGGCCGGCGGCGGCGCTGCTGCAAGCTT"
    "ACGTGGCATCACGATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACG"
    "TTCGGCCATGGACGTAGGCATCACGTGGCATCACGTTTAGGCCATGGACGTGGCATCACGA"
)

GC_RICH = ("GCCGGCGGCGGCGCTGCTGCAGCCGGCGCTGGCGGCTGCCGGCGGAGGCGCTGCTGCAGCC" * 6)


@pytest.mark.parametrize("template", [TEMPLATE, GC_RICH], ids=["mixed", "gc_rich"])
def test_no_returned_pair_forms_a_3_prime_dimer(client: TestClient, template: str) -> None:
    """The property the design path could not previously guarantee.

    GC-rich template included deliberately: that is where a 3' dimer becomes
    likely, and the issue notes the audit corpus happened to contain none, which
    reflected the corpus rather than the code.
    """
    r = client.post("/primers/design", json={
        "template": template, "tm_min": 45.0, "tm_max": 80.0,
        "gc_min": 0.0, "gc_max": 100.0,
        "product_min": 80, "product_max": 250, "max_pairs": 20,
    })
    assert r.status_code == 200, r.text
    pairs = r.json()
    if not pairs:
        pytest.skip("no pairs pass the other filters on this template")

    offenders = [
        (p["forward"]["sequence"], p["reverse"]["sequence"])
        for p in pairs
        if forms_3prime_dimer(p["forward"]["sequence"], p["reverse"]["sequence"])
    ]
    assert not offenders, f"{len(offenders)} of {len(pairs)} pairs form a 3' dimer: {offenders[:3]}"


def test_pairs_are_still_returned_for_an_ordinary_template(client: TestClient) -> None:
    """The check must not reject everything — that would be #62's hairpin again."""
    r = client.post("/primers/design", json={
        "template": TEMPLATE, "tm_min": 45.0, "tm_max": 75.0,
        "product_min": 80, "product_max": 250, "max_pairs": 20,
    })
    assert r.status_code == 200, r.text
    assert r.json(), "the dimer check rejected every pair"
