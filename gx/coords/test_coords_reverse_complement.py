"""CO-12 — reverse-complementing a record must not change what a feature reads.

Flipping a molecule end-for-end does not edit a gene: the same bases are still
in it, in the same order, read from the other strand. So the test is not "do
the numbers look mirrored" but "does the feature still extract the same
string". A single-part feature passes either way; the ones that pin the
convention are the two compound features in NC_001699.1.

Subject: ``backend/services/annotations.py:flip_annotations``, called from
``gui/panels/sequence_panel.py:_commit_edit(flip=True)`` for the Reverse
Complement operation.
"""

from __future__ import annotations

import pytest

from conftest import (
    bio_extract,
    feature_lines,
    flipped_sequence,
    requirement,
    reread,
    spliced,
)

CO_12 = requirement("CO-12")
LENGTH = CO_12["parameters"]["record_length"]
COMPOUND = CO_12["parameters"]["compound_features"]


@pytest.fixture
def flipped(jcv):
    return flipped_sequence(jcv)


def test_co_12_no_feature_is_lost_and_the_bases_are_the_reverse_complement(jcv, flipped):
    from Bio.Seq import Seq

    assert len(flipped.annotations) == CO_12["parameters"]["feature_count"]
    assert flipped.seq == str(Seq(jcv.seq).reverse_complement())
    assert flipped.length == LENGTH


def test_co_12_every_feature_still_extracts_the_same_string(jcv_bio, flipped):
    """The whole requirement in one assertion, re-sliced from the flipped
    string by index rather than trusting the app's own extraction."""
    want = [bio_extract(f, jcv_bio) for f in jcv_bio.features]
    wrong = []
    for i, ann in enumerate(flipped.annotations):
        got = spliced(flipped.seq, [tuple(p) for p in ann.parts], ann.strand.value)
        if got != want[i]:
            wrong.append((i, ann.feature_type, want[i][:24], got[:24]))
    assert wrong == [], f"{len(wrong)} of {len(want)} features read different bases: {wrong}"


def test_co_12_the_same_holds_after_export_and_re_read(jcv_bio, flipped):
    """Not just in memory: the written file has to say it too."""
    from backend.formats.genbank import write_genbank

    exported = write_genbank([flipped])
    after = reread(exported)
    wrong = [
        (i, f.type, str(f.location))
        for i, (f, g) in enumerate(zip(jcv_bio.features, after.features))
        if bio_extract(f, jcv_bio) != bio_extract(g, after)
    ]
    assert wrong == [], wrong


@pytest.mark.parametrize("case", COMPOUND, ids=lambda c: f"{c['feature_index']}-{c['type']}")
def test_co_12_compound_parts_stay_in_reading_order(flipped, case):
    """[s, e) -> [L-e, L-s) preserves the 5'-to-3' order of the part list, so
    the exon the feature reads first before the flip is still first after it."""
    ann = flipped.annotations[case["feature_index"]]
    assert ann.strand.value == case["strand_after"]
    assert [list(p) for p in ann.parts] == case["parts_after_expected"]


@pytest.mark.parametrize("case", COMPOUND, ids=lambda c: f"{c['feature_index']}-{c['type']}")
def test_co_12_compound_exports_with_its_exons_in_the_right_order(flipped, case):
    from backend.formats.genbank import write_genbank

    _key, loc = feature_lines(write_genbank([flipped]))[case["feature_index"]]
    assert loc == case["exported_after_expected"]


def test_co_12_spliced_cds_still_starts_at_its_start_codon(jcv_bio, flipped):
    """The concrete consequence: exons out of order means the CDS no longer
    begins where the protein does."""
    case = COMPOUND[0]
    ann = flipped.annotations[case["feature_index"]]
    got = spliced(flipped.seq, [tuple(p) for p in ann.parts], ann.strand.value)
    assert len(got) == case["extract_length"]
    assert got[:60] == case["extract_first_60"]


def test_co_12_origin_spanning_feature_keeps_its_bases(flipped):
    case = COMPOUND[1]
    ann = flipped.annotations[case["feature_index"]]
    got = spliced(flipped.seq, [tuple(p) for p in ann.parts], ann.strand.value)
    assert got == case["extract"]


@pytest.mark.parametrize("i", [0, 3, 5, 11, 16, 25, 26, 27])
def test_co_12_single_part_features_land_on_the_mirrored_interval(jcv, flipped, i):
    """[s, e) -> [L - e, L - s), with the strand inverted."""
    before = [tuple(p) for p in jcv.annotations[i].parts]
    assert len(before) == 1
    s, e = before[0]
    assert [tuple(p) for p in flipped.annotations[i].parts] == [(LENGTH - e, LENGTH - s)]
    flip_strand = {"+": "-", "-": "+", ".": "."}[jcv.annotations[i].strand.value]
    assert flipped.annotations[i].strand.value == flip_strand


def test_co_12_matches_biopythons_own_reverse_complement(jcv_bio, flipped):
    """The most direct oracle available: Biopython flips a record's features
    itself (``SeqRecord.reverse_complement(features=True)``). Compare the whole
    feature set, so neither list order nor a lucky index can hide a difference."""
    rc = jcv_bio.reverse_complement(features=True)
    want = sorted(
        (f.type, tuple(tuple(map(int, (p.start, p.end))) for p in f.location.parts),
         {1: "+", -1: "-"}.get(f.location.strand, "."))
        for f in rc.features
    )
    got = sorted(
        (a.feature_type, tuple(tuple(p) for p in a.parts), a.strand.value)
        for a in flipped.annotations
    )
    assert got == want


def test_co_12_flipping_twice_is_the_identity(jcv):
    """A control, and the reason a round-trip test cannot find an order error:
    reversing a part list twice restores it whether or not reversing was right."""
    twice = flipped_sequence(flipped_sequence(jcv))
    assert twice.seq == jcv.seq
    assert [([tuple(p) for p in a.parts], a.strand) for a in twice.annotations] == [
        ([tuple(p) for p in a.parts], a.strand) for a in jcv.annotations
    ]
