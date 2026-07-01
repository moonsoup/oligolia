"""Tests for off-target prediction (backend/crispr_offtarget.py)."""

from ..crispr_offtarget import (
    _mit_hit_score,
    _reverse_complement,
    scan_off_targets,
)


def test_reverse_complement() -> None:
    assert _reverse_complement("ATCG") == "CGAT"
    assert _reverse_complement("AAAA") == "TTTT"


def test_mit_hit_score_exact_is_one() -> None:
    assert _mit_hit_score([]) == 1.0


def test_mit_hit_score_decreases_with_mismatches() -> None:
    one = _mit_hit_score([1])
    two = _mit_hit_score([1, 2])
    three = _mit_hit_score([1, 2, 3])
    assert 0.0 < three < two < one <= 1.0


def test_seed_mismatch_penalized_more_than_distal() -> None:
    """A mismatch in the PAM-proximal seed hurts more than a distal one."""
    distal = _mit_hit_score([1])   # 5' end, tolerant
    seed = _mit_hit_score([20])    # PAM-proximal, intolerant
    assert seed < distal


def test_on_target_excluded() -> None:
    """A guide scanned against only its own origin has zero off-targets."""
    guide = "ACGTACGTACGTACGTACGT"
    reference = guide + "GG"  # protospacer + GG PAM, single on-target site
    result = scan_off_targets(guide, [reference], cas_family="cas9")
    assert result.summary["0"] == 0
    assert result.total == 0
    assert result.specificity_score == 100.0


def test_exact_off_target_detected() -> None:
    """A second identical protospacer+PAM is a real off-target and tanks score."""
    guide = "ACGTACGTACGTACGTACGT"
    site = guide + "GG"
    reference = site + "AAAAAAAA" + site  # on-target + one exact duplicate
    result = scan_off_targets(guide, [reference], cas_family="cas9")
    assert result.summary["0"] == 1
    assert result.total >= 1
    assert result.specificity_score < 100.0


def test_mismatch_off_target_bucketed() -> None:
    """A near-match with a single mismatch lands in the 1-mismatch bucket."""
    guide = "ACGTACGTACGTACGTACGT"
    on_target = guide + "GG"
    # Same protospacer with one substitution at position 5 (A->T), still + GG.
    mm_site = "ACGTTCGTACGTACGTACGT" + "GG"
    reference = on_target + "CCCCCC" + mm_site
    result = scan_off_targets(guide, [reference], cas_family="cas9")
    assert result.summary["0"] == 0  # on-target removed
    assert result.summary["1"] >= 1
    assert result.specificity_score < 100.0


def test_max_mismatches_bounds_buckets() -> None:
    guide = "ACGTACGTACGTACGTACGT"
    reference = guide + "GG"
    result = scan_off_targets(guide, [reference], cas_family="cas9", max_mismatches=2)
    assert set(result.summary.keys()) == {"0", "1", "2"}
