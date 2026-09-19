"""CR-6 — guide_length, honoured or defaulted, across every nuclease.

#64.3 was `guide_length` ignored by every branch. The oracle is the same
independent both-strand scan at the requested length, so this checks both that
the returned guides are the right size and that the search space moved with them.
"""

from __future__ import annotations

import pytest

from crispr_oracle import CANONICAL_LENGTH, requirement, scan

CAS_TYPES = ["SpCas9", "SpCas9-HF1", "AsCas12a", "LwaCas13a"]
LENGTHS = list(range(17, 25))


@pytest.mark.parametrize("cas", CAS_TYPES)
def test_cr_6_canonical_length_applies_when_the_caller_sets_none(puc19, design, cas):
    canonical = requirement("CR-6")["parameters"]["canonical_when_unset"]
    assert CANONICAL_LENGTH[cas] == canonical[cas]

    resp = design(puc19, cas, max_guides=50)
    assert {len(g.sequence) for g in resp.guides} == {canonical[cas]}
    assert resp.total_candidates == len(scan(cas, puc19, canonical[cas]))


@pytest.mark.parametrize("cas", CAS_TYPES)
@pytest.mark.parametrize("length", LENGTHS)
def test_cr_6_explicit_guide_length_is_honoured(puc19, design, cas, length):
    resp = design(puc19, cas, guide_length=length, max_guides=50)
    assert {len(g.sequence) for g in resp.guides} == {length}
    assert resp.total_candidates == len(scan(cas, puc19, length))


@pytest.mark.parametrize("cas", CAS_TYPES)
def test_cr_6_length_change_moves_the_search_at_real_scale(col1a1, design, cas):
    """The same check on the 24.5 kb human gene, where a stale length hides."""
    for length in (18, 21, 24):
        resp = design(col1a1, cas, guide_length=length, max_guides=50)
        assert {len(g.sequence) for g in resp.guides} == {length}
        assert resp.total_candidates == len(scan(cas, col1a1, length))
