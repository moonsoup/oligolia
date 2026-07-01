"""Tests for the Sequence model, including topology (is_circular)."""

from ..models.sequence import MoleculeType, Sequence


def test_is_circular_defaults_false() -> None:
    s = Sequence(id="x", seq="ATGC")
    assert s.is_circular is False


def test_is_circular_roundtrips_through_dump_and_validate() -> None:
    s = Sequence(id="p1", seq="ATGCATGC", molecule_type=MoleculeType.DNA, is_circular=True)
    dumped = s.model_dump()
    assert dumped["is_circular"] is True
    restored = Sequence.model_validate(dumped)
    assert restored.is_circular is True
    assert restored.seq == s.seq


def test_length_still_autofills_with_new_field() -> None:
    s = Sequence(id="x", seq="ATGCATGC", is_circular=True)
    assert s.length == 8  # model_post_init unaffected by the added field
