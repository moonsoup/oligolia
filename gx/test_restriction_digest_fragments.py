"""Digest fragments for linear and circular templates: RE-6, RE-7, RE-8.

Oracle: ``Bio.Restriction.<enz>.catalyse(seq, linear=...)`` — Biopython's own
fragmenting, which for a circular molecule bridges the last site back to the
first — plus fragment sequences re-derived by slicing the pinned template.
"""

from __future__ import annotations

import pytest
from Bio import Restriction
from Bio.Seq import Seq

from backend.routers.primers import _CURATED_ENZYMES, DigestRequest, digest

PANEL = list(_CURATED_ENZYMES)

# RE-6 / RE-7, recorded literals (register.json parameters.cases).
LINEAR_CASES = {
    "EcoRI": [396, 2290],
    "PvuII": [308, 322, 2056],
    "AvaII": [222, 627, 1837],
}
CIRCULAR_CASES = {
    "EcoRI": [2686],
    "PvuII": [322, 2364],
    "AvaII": [222, 2464],
}
PANEL_LINEAR_LENGTHS = [2, 2, 3, 6, 6, 6, 6, 10, 10, 88, 124, 183, 184, 222, 627, 1207]
PANEL_CIRCULAR_LENGTHS = [2, 2, 3, 6, 6, 6, 6, 10, 10, 88, 124, 183, 222, 811, 1207]

NO_SITE_ENZYMES = ["NotI", "XhoI", "NcoI", "ClaI", "EcoRV", "MluI", "NheI", "AgeI", "BglII", "MfeI"]


def _biopython_lengths(template: str, enzyme: str, *, linear: bool) -> list[int]:
    enz = getattr(Restriction, enzyme)
    return sorted(len(f) for f in enz.catalyse(Seq(template), linear=linear))


@pytest.mark.parametrize("enzyme", sorted(LINEAR_CASES))
def test_re_6_linear_fragments_match_biopython(puc19, enzyme):
    """RE-6: n cuts -> n+1 fragments, lengths equal catalyse(linear=True), sum = 2686."""
    expected = _biopython_lengths(puc19, enzyme, linear=True)
    assert expected == LINEAR_CASES[enzyme]
    n_cuts = len(getattr(Restriction, enzyme).search(Seq(puc19), linear=True))
    assert len(expected) == n_cuts + 1

    result = digest(DigestRequest(template=puc19, enzymes=[enzyme], is_circular=False))
    assert sorted(f.length for f in result.fragments) == expected
    assert sum(f.length for f in result.fragments) == len(puc19)
    assert result.template_length == len(puc19)
    for frag in result.fragments:
        assert frag.sequence == puc19[frag.start : frag.end]
        assert frag.length == len(frag.sequence)


def test_re_6_linear_whole_panel(puc19):
    """RE-6: the 15-cut panel digest of linear pUC19 gives 16 fragments summing to 2686."""
    expected_cuts = sorted(
        {p - 1 for name in PANEL for p in getattr(Restriction, name).search(Seq(puc19), linear=True)}
    )
    assert len(expected_cuts) == 15

    result = digest(DigestRequest(template=puc19, enzymes=PANEL, is_circular=False))
    assert len(result.fragments) == len(expected_cuts) + 1
    assert sorted(f.length for f in result.fragments) == PANEL_LINEAR_LENGTHS
    assert sum(f.length for f in result.fragments) == 2686

    # Reassembling the fragments in coordinate order must give the template back.
    in_order = sorted(result.fragments, key=lambda f: f.start)
    assert "".join(f.sequence for f in in_order) == puc19


@pytest.mark.parametrize("enzyme", sorted(CIRCULAR_CASES))
def test_re_7_circular_fragments_match_biopython(puc19, enzyme):
    """RE-7: n cuts -> n fragments (not n+1), lengths equal catalyse(linear=False)."""
    expected = _biopython_lengths(puc19, enzyme, linear=False)
    assert expected == CIRCULAR_CASES[enzyme]
    n_cuts = len(getattr(Restriction, enzyme).search(Seq(puc19), linear=False))
    assert len(expected) == n_cuts

    result = digest(DigestRequest(template=puc19, enzymes=[enzyme], is_circular=True))
    assert len(result.fragments) == n_cuts
    assert sorted(f.length for f in result.fragments) == expected
    assert sum(f.length for f in result.fragments) == len(puc19)
    for frag in result.fragments:
        expected_seq = (
            puc19[frag.start : frag.end]
            if frag.end > frag.start
            else puc19[frag.start :] + puc19[: frag.end]
        )
        assert frag.sequence == expected_seq
        assert frag.length == len(expected_seq)


def test_re_7_circular_whole_panel(puc19):
    """RE-7: the same 15 cuts give 15 fragments on the circle, one fewer than linear."""
    expected_cuts = sorted(
        {p - 1 for name in PANEL for p in getattr(Restriction, name).search(Seq(puc19), linear=False)}
    )
    assert len(expected_cuts) == 15

    result = digest(DigestRequest(template=puc19, enzymes=PANEL, is_circular=True))
    assert len(result.fragments) == 15
    assert sorted(f.length for f in result.fragments) == PANEL_CIRCULAR_LENGTHS
    assert sum(f.length for f in result.fragments) == 2686


@pytest.mark.parametrize("enzyme", NO_SITE_ENZYMES)
@pytest.mark.parametrize("circular", [False, True])
def test_re_8_enzyme_with_no_site(puc19, enzyme, circular):
    """RE-8: no site -> no cut, one fragment that is the whole template."""
    from backend.routers.primers import RestrictionRequest, restriction_sites

    assert getattr(Restriction, enzyme).search(Seq(puc19), linear=not circular) == []

    reported = {
        s.enzyme
        for s in restriction_sites(RestrictionRequest(template=puc19, is_circular=circular))
    }
    assert enzyme not in reported

    result = digest(DigestRequest(template=puc19, enzymes=[enzyme], is_circular=circular))
    assert result.cut_positions == []
    assert len(result.fragments) == 1
    assert result.fragments[0].length == 2686
    assert result.fragments[0].sequence == puc19
