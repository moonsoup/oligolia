"""A recognition site that straddles the origin of a circular template: RE-5.

pUC19 rotated left by 398 bases puts the unique EcoRI site across the junction:
the rotated template ends with GAA and begins with TTC. Oracle: the doubled-
sequence window search re-derived by index, and ``EcoRI.search(seq, linear=False)``.
"""

from __future__ import annotations

from Bio import Restriction
from Bio.Seq import Seq
from conftest import find_starts, find_starts_circular

from backend.routers.primers import (
    DigestRequest,
    RestrictionRequest,
    digest,
    restriction_sites,
)

ROTATION = 398


def _rotated(puc19: str) -> str:
    return puc19[ROTATION:] + puc19[:ROTATION]


def test_re_5_origin_spanning_site_found_exactly_once(puc19):
    """RE-5: the straddling site is reported once, at its true start 2683."""
    rot = _rotated(puc19)
    assert len(rot) == 2686
    assert rot[-3:] + rot[:3] == "GAATTC"

    # Independent oracles: index-derived circular windows, and Biopython.
    assert find_starts_circular(rot, Restriction.EcoRI.site) == [2683]
    assert find_starts(rot, Restriction.EcoRI.site) == []
    assert Restriction.EcoRI.search(Seq(rot), linear=False) == [2685]

    sites = restriction_sites(
        RestrictionRequest(template=rot, enzymes=["EcoRI"], is_circular=True)
    )
    assert len(sites) == 1
    assert sites[0].positions == [2683]
    assert sites[0].count == 1


def test_re_5_origin_spanning_site_absent_when_linear(puc19):
    """RE-5: the same template read as linear has no EcoRI site and one fragment."""
    rot = _rotated(puc19)
    assert Restriction.EcoRI.search(Seq(rot), linear=True) == []

    sites = restriction_sites(
        RestrictionRequest(template=rot, enzymes=["EcoRI"], is_circular=False)
    )
    assert sites == []

    result = digest(DigestRequest(template=rot, enzymes=["EcoRI"], is_circular=False))
    assert result.cut_positions == []
    assert [f.length for f in result.fragments] == [2686]
    assert result.fragments[0].sequence == rot


def test_re_5_origin_spanning_cut_and_single_fragment(puc19):
    """RE-5: one cut at 2684 (Bio 2685 - 1), one fragment of the full 2686 bp."""
    rot = _rotated(puc19)
    expected_cut = Restriction.EcoRI.search(Seq(rot), linear=False)[0] - 1
    assert expected_cut == 2684
    assert 2683 + Restriction.EcoRI.fst5 == expected_cut

    result = digest(DigestRequest(template=rot, enzymes=["EcoRI"], is_circular=True))
    assert result.cut_positions == [expected_cut]
    assert len(result.fragments) == 1
    frag = result.fragments[0]
    assert frag.length == 2686
    assert frag.sequence == rot[expected_cut:] + rot[:expected_cut]
    # A circular one-cut digest linearises the plasmid: the fragment carries the
    # same bases as the template, read from the cut.
    assert sorted(frag.sequence) == sorted(rot)


def test_re_5_origin_spanning_site_found_once_for_every_rotation(puc19):
    """RE-5: for every rotation the circular search finds the site exactly once."""
    n = len(puc19)
    for rotation in (0, 1, 393, 396, 397, 398, 399, 400, 2685):
        rot = puc19[rotation:] + puc19[:rotation]
        expected = find_starts_circular(rot, Restriction.EcoRI.site)
        assert expected == [(395 - rotation) % n]
        sites = restriction_sites(
            RestrictionRequest(template=rot, enzymes=["EcoRI"], is_circular=True)
        )
        assert [s.positions for s in sites] == [expected], rotation
