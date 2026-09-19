"""Cut-position convention, overlapping sites, coinciding cuts: RE-4, RE-9, RE-14.

Oracles: ``Bio.Restriction.<enz>.search`` (1-based, "first base of the 3'
fragment"), the ``fst5`` offset from the REBASE table, and positions re-derived
by index from the pinned sequence.
"""

from __future__ import annotations

import pytest
from Bio import Restriction
from Bio.Seq import Seq
from conftest import find_starts

from backend.routers.primers import (
    DigestRequest,
    RestrictionRequest,
    digest,
    restriction_sites,
)

# RE-4, recorded literals (register.json parameters.cases).
CUT_CASES = {
    "EcoRI": {"site_start": 395, "fst5": 1, "cut": 396, "bio_1based": 397},
    "NdeI": {"site_start": 182, "fst5": 2, "cut": 184, "bio_1based": 185},
    "SmaI": {"site_start": 411, "fst5": 3, "cut": 414, "bio_1based": 415},
    "KpnI": {"site_start": 407, "fst5": 5, "cut": 412, "bio_1based": 413},
    "PvuII": {"site_start": 305, "fst5": 3, "cut": 308, "bio_1based": 309},
}


@pytest.mark.parametrize("enzyme", sorted(CUT_CASES))
def test_re_4_cut_position_convention(puc19, enzyme):
    """RE-4: cut = Bio 1-based position - 1 = site start + fst5, 0-based, top strand."""
    case = CUT_CASES[enzyme]
    enz = getattr(Restriction, enzyme)

    assert enz.fst5 == case["fst5"]
    bio = enz.search(Seq(puc19), linear=True)
    assert bio[0] == case["bio_1based"]
    expected_cut = case["bio_1based"] - 1
    assert expected_cut == case["cut"]

    site_start = find_starts(puc19, enz.site)[0]
    assert site_start == case["site_start"]
    assert site_start + enz.fst5 == expected_cut

    # The site really is where the oracle says it is, and the cut falls fst5
    # bases into it.
    assert puc19[site_start : site_start + enz.size] == enz.site

    result = digest(DigestRequest(template=puc19, enzymes=[enzyme]))
    assert result.cut_positions[0] == expected_cut
    # "first base of the 3' fragment": exactly one fragment starts at the cut,
    # and it is the template slice beginning there.
    downstream = [f for f in result.fragments if f.start == expected_cut]
    assert len(downstream) == 1
    frag = downstream[0]
    assert frag.sequence == puc19[frag.start : frag.end]


def test_re_4_whole_panel_cut_positions(puc19):
    """RE-4: the convention holds for every cutting panel enzyme at once."""
    from backend.routers.primers import _CURATED_ENZYMES

    expected = sorted(
        {p - 1 for name in _CURATED_ENZYMES for p in getattr(Restriction, name).search(Seq(puc19))}
    )
    assert expected == [
        184, 308, 396, 406, 412, 414, 417, 423, 429, 439, 445, 447, 630, 1837, 2059,
    ]
    got = digest(DigestRequest(template=puc19, enzymes=list(_CURATED_ENZYMES)))
    assert got.cut_positions == expected


def test_re_9_overlapping_sites_are_each_reported():
    """RE-9: two ClaI sites sharing two bases are both found."""
    template = "AAATCGATCGATTTT"
    assert find_starts(template, Restriction.ClaI.site) == [2, 6]
    bio = Restriction.ClaI.search(Seq(template), linear=True)
    assert bio == [5, 9]
    expected_cuts = [p - 1 for p in bio]
    assert expected_cuts == [4, 8]

    sites = restriction_sites(RestrictionRequest(template=template, enzymes=["ClaI"]))
    assert [s.positions for s in sites] == [[2, 6]]

    result = digest(DigestRequest(template=template, enzymes=["ClaI"]))
    assert result.cut_positions == expected_cuts
    assert sorted(f.length for f in result.fragments) == [4, 4, 7]
    assert sum(f.length for f in result.fragments) == len(template)


def test_re_9_overlapping_sites_of_different_enzymes(puc19):
    """RE-9: KpnI 407-412, SmaI 411-416 and AvaI 411-416 overlap in the pUC19 MCS."""
    assert puc19[407:413] == "GGTACC"
    assert puc19[411:417] == "CCCGGG"

    sites = {
        s.enzyme: list(s.positions)
        for s in restriction_sites(
            RestrictionRequest(template=puc19, enzymes=["KpnI", "SmaI", "AvaI"])
        )
    }
    assert sites == {"KpnI": [407], "SmaI": [411], "AvaI": [411]}

    expected_cuts = sorted(
        {p - 1 for e in ("KpnI", "SmaI", "AvaI") for p in getattr(Restriction, e).search(Seq(puc19))}
    )
    assert expected_cuts == [412, 414]
    result = digest(DigestRequest(template=puc19, enzymes=["KpnI", "SmaI", "AvaI"]))
    assert result.cut_positions == expected_cuts


def test_re_14_digest_is_independent_of_enzyme_order(puc19):
    """RE-14: AvaI and KpnI both cut at 412 with opposite overhang chemistry."""
    assert Restriction.AvaI.ovhg == -4
    assert Restriction.KpnI.ovhg == 4
    assert Restriction.AvaI.search(Seq(puc19))[0] - 1 == 412
    assert Restriction.KpnI.search(Seq(puc19))[0] - 1 == 412

    one = digest(DigestRequest(template=puc19, enzymes=["AvaI", "KpnI"]))
    two = digest(DigestRequest(template=puc19, enzymes=["KpnI", "AvaI"]))

    assert [f.model_dump() for f in one.fragments] == [f.model_dump() for f in two.fragments]
