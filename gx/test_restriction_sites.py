"""Recognition-site finding: RE-1, RE-2, RE-3, RE-10.

Oracles: occurrences re-derived by index from the pinned pUC19 bytes, the
REBASE recognition/cut data read straight off ``Bio.Restriction`` (pinned as
gx/corpus/biopython-1.85/Bio/Restriction/Restriction_Dictionary.py), and
``Bio.Restriction.<enz>.search`` itself. Never the subject's own output.
"""

from __future__ import annotations

from Bio import Restriction
from Bio.Seq import Seq
from conftest import (
    find_starts,
    iupac_regex,
    iupac_reverse_complement,
)

from backend.routers.primers import (
    _CURATED_ENZYMES,
    RestrictionRequest,
    restriction_sites,
)

# RE-1, recorded literals (register.json parameters.expected_site_starts).
PUC19_SITE_STARTS = {
    "NdeI": [182],
    "PvuII": [305, 627],
    "EcoRI": [395],
    "SacI": [401],
    "KpnI": [407],
    "SmaI": [411],
    "AvaI": [411],
    "BamHI": [416],
    "XbaI": [422],
    "SalI": [428],
    "PstI": [434],
    "SphI": [440],
    "HindIII": [446],
    "AvaII": [1836, 2058],
}


def _app_sites(template: str, *, circular: bool = False, enzymes=None) -> dict[str, list[int]]:
    result = restriction_sites(
        RestrictionRequest(template=template, enzymes=enzymes, is_circular=circular)
    )
    return {s.enzyme: list(s.positions) for s in result}


def test_re_1_site_starts_on_linear_puc19(puc19):
    """RE-1: reported positions are the 0-based starts of every real occurrence."""
    oracle = {}
    for name in _CURATED_ENZYMES:
        starts = find_starts(puc19, getattr(Restriction, name).site)
        if starts:
            oracle[name] = starts

    # The index-derived oracle must itself match the values recorded in the
    # register, so a drift in the pinned input cannot silently rewrite the test.
    assert oracle == PUC19_SITE_STARTS

    assert _app_sites(puc19) == oracle


def test_re_1_count_and_pattern_fields(puc19):
    """RE-1: ``count`` is the number of occurrences, ``cut_pattern`` the REBASE site."""
    for site in restriction_sites(RestrictionRequest(template=puc19)):
        enz = getattr(Restriction, site.enzyme)
        assert site.cut_pattern == enz.site
        assert site.count == len(find_starts(puc19, enz.site))


def test_re_2_panel_recognition_sites_are_palindromic():
    """RE-2: every panel site equals its own IUPAC reverse complement."""
    not_palindromic = {
        name: getattr(Restriction, name).site
        for name in _CURATED_ENZYMES
        if getattr(Restriction, name).site
        != iupac_reverse_complement(getattr(Restriction, name).site)
    }
    assert not_palindromic == {}


def test_re_2_reverse_strand_gives_mirror_positions(puc19):
    """RE-2: searching the complementary strand yields n-(start+size)."""
    rc = str(Seq(puc19).reverse_complement())
    forward = _app_sites(puc19)
    reverse = _app_sites(rc)
    n = len(puc19)

    expected = {
        name: sorted(n - (p + getattr(Restriction, name).size) for p in starts)
        for name, starts in PUC19_SITE_STARTS.items()
    }
    assert forward == PUC19_SITE_STARTS
    assert reverse == expected


def test_re_3_site_on_bottom_strand_is_found(puc19, panel_with):
    """RE-3: BsaI's site occurs in pUC19 only on the bottom strand."""
    panel_with("BsaI")
    enz = Restriction.BsaI

    top = find_starts(puc19, enz.site)
    bottom = find_starts(puc19, iupac_reverse_complement(enz.site))
    assert top == []            # GGTCTC does not occur on the top strand
    assert bottom == [1765]     # GAGACC does, so the site is on the bottom strand
    assert len(enz.search(Seq(puc19), linear=True)) == 1

    assert _app_sites(puc19, enzymes=["BsaI"]) == {"BsaI": [1765]}


def test_re_3_site_count_matches_cut_count(puc19, panel_with):
    """RE-3: the site finder and the digest must agree on how often BsaI cuts."""
    from backend.routers.primers import DigestRequest, digest

    panel_with("BsaI")
    cuts = digest(DigestRequest(template=puc19, enzymes=["BsaI"])).cut_positions
    assert cuts == [1760]  # Bio.Restriction: 1-based 1761 -> 0-based 1760

    sites = _app_sites(puc19, enzymes=["BsaI"]).get("BsaI", [])
    assert len(sites) == len(cuts)


def test_re_10_ambiguous_avai_all_four_resolutions():
    """RE-10: CYCGRG matches CCCGAG, CTCGAG, CCCGGG and CTCGGG."""
    template = "AAA" + "CCCGAG" + "TT" + "CTCGAG" + "TT" + "CCCGGG" + "TT" + "CTCGGG" + "AAA"
    assert iupac_regex(Restriction.AvaI.site) == "C[CT]CG[AG]G"
    assert find_starts(template, Restriction.AvaI.site) == [3, 11, 19, 27]
    assert [p - 1 - Restriction.AvaI.fst5 for p in Restriction.AvaI.search(Seq(template))] == [
        3,
        11,
        19,
        27,
    ]

    assert _app_sites(template, enzymes=["AvaI"]) == {"AvaI": [3, 11, 19, 27]}


def test_re_10_ambiguous_avaii_both_resolutions():
    """RE-10: GGWCC matches GGACC and GGTCC."""
    template = "TT" + "GGACC" + "TTT" + "GGTCC" + "TT"
    assert iupac_regex(Restriction.AvaII.site) == "GG[AT]CC"
    assert find_starts(template, Restriction.AvaII.site) == [2, 10]

    assert _app_sites(template, enzymes=["AvaII"]) == {"AvaII": [2, 10]}


def test_re_10_ambiguous_enzymes_on_puc19(puc19):
    """RE-10: the degenerate panel enzymes find their real pUC19 occurrences."""
    assert find_starts(puc19, Restriction.AvaI.site) == [411]
    assert find_starts(puc19, Restriction.AvaII.site) == [1836, 2058]

    got = _app_sites(puc19, enzymes=["AvaI", "AvaII"])
    assert got == {"AvaI": [411], "AvaII": [1836, 2058]}
