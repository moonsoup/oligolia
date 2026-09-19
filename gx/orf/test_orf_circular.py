"""OR-12 — a circle has no beginning, so the ORF set cannot depend on one.

pUC19 (GenBank L09137.2, 2686 bp, circular) carries the bla ampicillin-resistance
gene as a complete 286-residue ORF at top-strand ``[1625:2486)`` in frame -3.
Where the flat file happens to start is an editorial choice about a circular
molecule; rotating it must not change what genes the molecule has.

No second ORF finder is needed. A rotation is the identity on a circle, so the
subject's own answer for the unrotated molecule is the expected value, and OR-7
has already pinned that answer against an independent finder.
"""

from __future__ import annotations

import pytest

ROTATIONS = [1700, 2000, 2400]
BLA_SPAN = (1625, 2486)
BLA_LENGTH_AA = 286
BLA_FRAME = -3


def test_or_12_bla_is_a_complete_orf_in_the_unrotated_record(puc19, find_orfs) -> None:
    """The premise, stated separately so a failure downstream is unambiguous."""
    orfs = find_orfs(puc19, min_length_aa=30).orfs
    bla = [o for o in orfs if (o.start, o.end) == BLA_SPAN]
    assert len(bla) == 1, f"no ORF at {BLA_SPAN} in pUC19"
    assert (bla[0].frame, bla[0].length_aa) == (BLA_FRAME, BLA_LENGTH_AA)


@pytest.mark.parametrize("offset", ROTATIONS)
def test_or_12_rotating_the_origin_does_not_lose_the_bla_gene(
    offset, puc19, find_orfs
) -> None:
    """Each rotation puts the origin inside bla's 1625..2486 span."""
    assert BLA_SPAN[0] < offset < BLA_SPAN[1]
    rotated = puc19[offset:] + puc19[:offset]
    assert len(rotated) == 2686

    lengths = sorted((o.length_aa for o in find_orfs(rotated, min_length_aa=30).orfs), reverse=True)
    assert BLA_LENGTH_AA in lengths, (
        f"pUC19 read from base {offset} reports no 286-residue ORF; the bla gene "
        f"now spans the origin. Longest ORF found: {lengths[0] if lengths else None}"
    )


@pytest.mark.parametrize("offset", ROTATIONS)
def test_or_12_the_orf_set_is_rotation_invariant(offset, puc19, find_orfs) -> None:
    """The multiset of proteins is the same molecule's, whatever base it starts at."""
    from collections import Counter

    rotated = puc19[offset:] + puc19[:offset]
    reference = Counter(o.protein for o in find_orfs(puc19, min_length_aa=30).orfs)
    result = Counter(o.protein for o in find_orfs(rotated, min_length_aa=30).orfs)

    lost = reference - result
    gained = result - reference
    assert (lost, gained) == (Counter(), Counter()), (
        f"rotating pUC19 by {offset}: {sum(lost.values())} ORF(s) lost "
        f"(longest {max((len(p) for p in lost), default=0)} aa), "
        f"{sum(gained.values())} gained "
        f"(longest {max((len(p) for p in gained), default=0)} aa)"
    )
