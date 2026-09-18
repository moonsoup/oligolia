"""Differential oracle for primer Tm — the first of the #75 oracle layer.

#50: `_tm_nearest_neighbor` was named for nearest-neighbour thermodynamics and
implemented the 1989 Wallace/GC-count rule. An independent audit over 120 primers
from 8 real hEDS-panel templates found 97 of them (81%) more than 3 degC out, worst
case +9.0 degC, systematically worse on GC-rich sequence — exactly where stacking
energy matters most.

The oracle is Biopython's `Bio.SeqUtils.MeltingTemp.Tm_NN`, already a dependency.

THE PARAMETER SET IS PART OF THE ASSERTION. Tm_NN takes four literature
nearest-neighbour tables, both strand concentrations, four salt species and one of
seven salt corrections; "Tm within 1 degC of Biopython" does not identify a number
without them. These are Biopython's documented defaults, which is also what the
audit used:

    nn_table = DNA_NN3 (Allawi & SantaLucia 1997, the "unified" set)
    saltcorr = 5 (Owczarzy et al. 2004)
    dnac1 = dnac2 = 25 nM,  Na = 50 mM,  K = Tris = Mg = dNTPs = 0

If those ever change, this file changes with them, deliberately and visibly.
"""

from __future__ import annotations

import pytest
from Bio.SeqUtils import MeltingTemp as mt

from backend.routers.primers import _tm_nearest_neighbor, _tm_wallace

TOLERANCE_C = 0.5

#: The pinned parameter set. Passed explicitly so the test does not silently
#: follow a future change to Biopython's defaults.
NN_PARAMS = dict(
    nn_table=mt.DNA_NN3,
    saltcorr=5,
    dnac1=25,
    dnac2=25,
    Na=50,
    K=0,
    Tris=0,
    Mg=0,
    dNTPs=0,
)


def oracle(seq: str) -> float:
    return mt.Tm_NN(seq, **NN_PARAMS)


# A GC sweep, all 20-mers, roughly 30% -> 80% GC.
SWEEP = [
    "ATATTATAGCATATTATAGC",  # very AT-rich
    "ATCAGTATCAGTATCAGTAT",
    "ACGTACGTACGTACGTACGT",  # 50%
    "ACGTGGCATCACGATGGCCT",
    "GCCTGTGGGCATTTGGCCAA",
    "GGCCGGCCGGCCGGCCAATT",
    "GGCCGGCCGGCCGGCCGGCC",  # 100% GC
]


@pytest.mark.parametrize("seq", SWEEP)
def test_tm_matches_biopython_nn_across_a_gc_sweep(seq: str) -> None:
    assert abs(_tm_nearest_neighbor(seq) - oracle(seq)) <= TOLERANCE_C, (
        seq, _tm_nearest_neighbor(seq), oracle(seq)
    )


def test_tm_sees_base_order_not_just_gc_count() -> None:
    """The Wallace rule's blind spot, stated as a property.

    These two 20-mers have identical length and identical GC count, so the old
    implementation returned the same Tm for both. Real nearest-neighbour
    thermodynamics cannot: adjacent-base stacking differs.
    """
    a = "GGGGGGGGGGAAAAAAAAAA"
    b = "GAGAGAGAGAGAGAGAGAGA"
    assert len(a) == len(b)
    assert a.count("G") + a.count("C") == b.count("G") + b.count("C")

    assert abs(_tm_nearest_neighbor(a) - _tm_nearest_neighbor(b)) > 1.0, (
        "same length and GC count gave the same Tm — this is the Wallace rule (#50)"
    )


def test_a_gc_rich_primer_is_not_under_predicted() -> None:
    """The shape of the worst case: claimed 57.6 degC, actual 66.6 degC."""
    gc_rich = "GCCGGCGGCGGCGCTGCTGC"
    assert abs(_tm_nearest_neighbor(gc_rich) - oracle(gc_rich)) <= TOLERANCE_C
    # And the old rule really was far out on it, so the test has teeth.
    assert _tm_wallace(gc_rich) < oracle(gc_rich) - 3.0, (
        "expected the Wallace rule to under-predict a GC-rich primer by >3 degC"
    )


def test_the_wallace_rule_is_still_available_and_unchanged() -> None:
    """Kept as an option, not deleted — just no longer misnamed.

    Values pinned from the original implementation so a future edit to
    _tm_wallace is a visible change rather than a silent one.
    """
    # len < 14 branch: 2*AT + 4*GC
    assert _tm_wallace("ACGTACGTACGT") == pytest.approx(2 * 6 + 4 * 6)
    # len >= 14 branch: 64.9 + 41 * (gc - 16.4) / (at + gc)
    seq = "ACGTACGTACGTACGTACGT"
    expected = 64.9 + 41 * (10 - 16.4) / 20
    assert _tm_wallace(seq) == pytest.approx(expected)


def test_tm_is_monotonic_in_gc_for_a_fixed_backbone() -> None:
    """A sanity property no reasonable Tm function should violate."""
    seqs = [
        "ATATATATATATATATATAT",
        "ATATATATATGCATATATAT",
        "ATATGCATATGCATATGCAT",
        "GCATGCATGCATGCATGCAT",
        "GCGCGCATGCGCGCATGCGC",
    ]
    tms = [_tm_nearest_neighbor(s) for s in seqs]
    assert tms == sorted(tms), tms


# --- end to end: what the endpoint reports, not just what the helper computes ---

def test_every_primer_the_endpoint_returns_has_an_nn_tm(client) -> None:
    """The check that would have caught #50 at the API level.

    The existing primer tests assert `tm_min <= tm <= tm_max`, which the filter
    guarantees by construction — they would pass for any Tm function at all,
    including the one that was 9 degC out. This compares the reported number
    against the oracle for the primer's own sequence.
    """
    template = (
        "ATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACGTGGCATCACGTTT"
        "CGGCCATGGACGTAGGCATCACGTGGCATCACGTGGCCGGCGGCGGCGCTGCTGCAAGCTT"
        "ACGTGGCATCACGATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACG"
    )
    r = client.post("/primers/design", json={
        "template": template, "tm_min": 45.0, "tm_max": 75.0,
        "product_min": 80, "product_max": 200, "max_pairs": 10,
    })
    assert r.status_code == 200, r.text
    pairs = r.json()
    assert pairs, "no primer pairs returned for this template"

    for pair in pairs:
        for role in ("forward", "reverse"):
            p = pair[role]
            expected = oracle(p["sequence"])
            assert abs(p["tm"] - expected) <= TOLERANCE_C, (
                role, p["sequence"], p["tm"], expected
            )


def test_the_endpoint_no_longer_reports_the_wallace_number(client) -> None:
    """Teeth: where the two rules disagree, the reported value must be the NN one.

    Uses the same template as the test above, which is known to yield pairs, so
    this cannot quietly turn into a skip — a test that stops running is how #50
    survived a green suite in the first place.
    """
    template = (
        "ATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACGTGGCATCACGTTT"
        "CGGCCATGGACGTAGGCATCACGTGGCATCACGTGGCCGGCGGCGGCGCTGCTGCAAGCTT"
        "ACGTGGCATCACGATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACG"
    )
    r = client.post("/primers/design", json={
        "template": template, "tm_min": 45.0, "tm_max": 75.0,
        "product_min": 80, "product_max": 200, "max_pairs": 10,
    })
    assert r.status_code == 200, r.text
    pairs = r.json()
    assert pairs, "no primer pairs returned for this template"

    checked = 0
    for pair in pairs:
        for role in ("forward", "reverse"):
            p_ = pair[role]
            wallace = _tm_wallace(p_["sequence"])
            nn = oracle(p_["sequence"])
            assert abs(p_["tm"] - nn) <= TOLERANCE_C, (p_["sequence"], p_["tm"], nn)
            if abs(wallace - nn) > 2.0:
                assert abs(p_["tm"] - wallace) > 1.0, (
                    f"reported Tm still matches the Wallace rule for {p_['sequence']}"
                )
                checked += 1
    assert checked, "no primer here distinguishes the two rules, so this proves nothing"
