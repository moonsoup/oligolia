"""Tests for the bundled reference feature library (issue #42)."""

from Bio.Seq import Seq

from ..formats import load_common_features


def test_library_loads_enough_parts() -> None:
    parts = load_common_features()
    assert len(parts) >= 15


def test_every_part_is_well_formed() -> None:
    parts = load_common_features()
    for p in parts:
        assert p.name
        assert p.feature_type
        assert p.sequence and set(p.sequence.upper()) <= set("ACGTUN")
        assert p.source, f"{p.name} must cite a source"


def test_part_names_are_unique() -> None:
    parts = load_common_features()
    names = [p.name for p in parts]
    assert len(names) == len(set(names))


def test_coding_parts_translate_to_documented_peptide() -> None:
    """Each coding part's DNA must actually encode the peptide it claims to."""
    parts = load_common_features()
    coding = [p for p in parts if p.translation]
    assert coding, "expected some coding parts (tags/cleavage sites)"
    for p in coding:
        assert len(p.sequence) % 3 == 0, f"{p.name}: not a whole number of codons"
        translated = str(Seq(p.sequence).translate())
        assert translated == p.translation, f"{p.name}: {translated} != {p.translation}"


def test_returns_independent_lists() -> None:
    """Mutating a returned list must not corrupt the cached source."""
    a = load_common_features()
    n = len(a)
    a.pop()
    assert len(load_common_features()) == n


def test_load_from_explicit_path(tmp_path) -> None:
    import json
    custom = tmp_path / "parts.json"
    custom.write_text(json.dumps({"parts": [
        {"name": "x", "feature_type": "misc", "sequence": "ATGC", "source": "test"}
    ]}))
    parts = load_common_features(custom)
    assert len(parts) == 1 and parts[0].name == "x"
