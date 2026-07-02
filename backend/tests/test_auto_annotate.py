"""Tests for homology-based auto-annotation (issue #43)."""

from Bio.Seq import Seq

from ..formats import auto_annotate, load_common_features

LIB = {p.name: p.sequence for p in load_common_features()}
FILLER = "GCTAGCTTACGATCGGATCACTGATCCAGTACGATCGATCACGTAGCATCGATCGTACG"


def _labels(anns) -> set[str]:
    return {a.qualifiers["label"] for a in anns}


def test_detects_exact_forward_parts() -> None:
    seq = FILLER + LIB["T7 promoter"] + FILLER + LIB["Myc tag"] + FILLER
    anns = auto_annotate(seq)
    assert {"T7 promoter", "Myc tag"} <= _labels(anns)
    for a in anns:
        assert a.qualifiers["auto_detected"] == "true"


def test_coordinates_point_at_the_match() -> None:
    seq = FILLER + LIB["T7 promoter"] + FILLER
    ann = next(a for a in auto_annotate(seq) if a.qualifiers["label"] == "T7 promoter")
    assert seq[ann.start:ann.end] == LIB["T7 promoter"]


def test_detects_reverse_strand_match() -> None:
    ha = LIB["HA tag"]
    seq = FILLER + str(Seq(ha).reverse_complement()) + FILLER
    hits = [a for a in auto_annotate(seq) if a.qualifiers["label"] == "HA tag"]
    assert hits, "HA tag on reverse strand not found"
    hit = hits[0]
    assert hit.strand.value == "-"
    # Forward-oriented coordinates re-derive the reverse-complemented match.
    assert str(Seq(seq[hit.start:hit.end]).reverse_complement()) == ha


def test_tolerates_a_mismatch() -> None:
    strep = LIB["Strep-tag II"]
    variant = "A" + strep[1:] if strep[0] != "A" else "C" + strep[1:]  # 1 substitution
    seq = FILLER + variant + FILLER
    hits = [a for a in auto_annotate(seq) if a.qualifiers["label"] == "Strep-tag II"]
    assert hits
    assert hits[0].qualifiers["mismatches"] == "1"


def test_rejects_too_many_mismatches() -> None:
    myc = LIB["Myc tag"]
    # Corrupt ~1/3 of positions — well beyond the 10% identity threshold.
    corrupted = "".join("A" if i % 3 == 0 else c for i, c in enumerate(myc))
    seq = FILLER + corrupted + FILLER
    assert "Myc tag" not in _labels(auto_annotate(seq))


def test_short_parts_below_min_length_are_skipped() -> None:
    # The 6 bp RBS should never be reported (below min_length), even if present.
    seq = FILLER + LIB["RBS (T7 gene 10)"] + FILLER
    assert "RBS (T7 gene 10)" not in _labels(auto_annotate(seq))


def test_clean_random_sequence_yields_no_hits() -> None:
    # Filler alone contains no library parts.
    assert auto_annotate(FILLER * 3) == []


def test_auto_annotate_finds_ampr_and_ori_in_puc19() -> None:
    """Acceptance (#43): scanning real pUC19 identifies AmpR and its ColE1 ori."""
    from ._puc19 import PUC19
    labels = {a.qualifiers["label"] for a in auto_annotate(PUC19)}
    assert "AmpR (bla)" in labels
    assert "ColE1 ori" in labels
