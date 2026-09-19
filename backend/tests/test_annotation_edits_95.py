"""#95: an edit must not leave a feature asserting something the bases deny.

Two ways `shift_annotations`' "contains the edit" branch used to keep a feature
it had no business keeping:

* deleting exactly a feature's own bases collapsed it to zero length, which the
  exporter wrote as the INSDC between-position `n^n+1` — a site *between* two
  bases — or, at the head of a record, as the impossible `0^1`;
* any edit inside a CDS left the pre-edit `/translation` in place, so the
  exported record declared a protein its own bases no longer encode. An
  equal-length replacement is the sharpest form: no coordinate moves at all.

The arithmetic for the cases that were already right lives in
`test_annotation_edits.py`; this file is only the two wrong turns and the guard
that the right turns stayed right.
"""

from __future__ import annotations

from Bio.Seq import Seq

from backend.formats.genbank import write_genbank
from backend.models.sequence import Annotation, MoleculeType, Sequence, Strand
from backend.services.annotations import shift_annotations

# M, ten alanines, stop — a complete little CDS sitting at [5, 41).
CDS_BASES = "ATG" + "GCT" * 10 + "TAA"
PROTEIN = "M" + "A" * 10
SEQ = "CCCCC" + CDS_BASES + "GGGGG"
CDS_START, CDS_END = 5, 5 + len(CDS_BASES)
SECOND_CODON = (CDS_START + 3, CDS_START + 6)  # the first GCT


def _ann(start: int, end: int, name: str = "f", parts=None, strand=Strand.PLUS,
         feature_type: str = "gene", **qualifiers) -> Annotation:
    return Annotation(
        feature_type=feature_type, start=start, end=end,
        parts=parts or [(start, end)], strand=strand,
        qualifiers={"gene": name, **qualifiers},
    )


def _cds(strand=Strand.PLUS, parts=None, protein: str = PROTEIN) -> Annotation:
    return _ann(CDS_START, CDS_END, "orf", parts=parts, strand=strand,
                feature_type="CDS", translation=protein, codon_start="1")


def _record(annotations: list[Annotation], seq: str = SEQ) -> Sequence:
    return Sequence(id="s", name="s", seq=seq, length=len(seq),
                    molecule_type=MoleculeType.DNA, annotations=annotations)


# ── nothing left to describe ────────────────────────────────────────────────

def test_deleting_exactly_a_features_bases_reports_it_lost() -> None:
    """It used to satisfy `p_start <= start and p_end >= end` and survive."""
    doomed = _ann(10, 20, "doomed")
    kept, dropped = shift_annotations([doomed], start=10, end=20, inserted=0)
    assert kept == []
    assert [a.qualifiers["gene"] for a in dropped] == ["doomed"]


def test_a_wholly_deleted_feature_is_not_exported_as_a_between_position() -> None:
    kept, _ = shift_annotations(
        [_ann(10, 20, "doomed"), _ann(30, 40, "after")], start=10, end=20, inserted=0
    )
    text = write_genbank([_record(kept, SEQ[:10] + SEQ[20:])])
    assert "^" not in text.split("ORIGIN")[0], text


def test_deleting_the_first_bases_does_not_export_a_location_at_base_zero() -> None:
    """A GenBank location counts from 1, so `0^1` is not a location at all."""
    kept, dropped = shift_annotations([_ann(0, 12, "head")], start=0, end=12, inserted=0)
    assert [a.qualifiers["gene"] for a in dropped] == ["head"]
    assert "0^1" not in write_genbank([_record(kept, SEQ[12:])])


def test_a_spliced_feature_whose_part_is_wholly_deleted_is_dropped() -> None:
    """One exon deleted is the same loss as half of one: the feature goes."""
    spliced = _ann(4, 40, "cds", parts=[(4, 10), (29, 40)])
    kept, dropped = shift_annotations([spliced], start=29, end=40, inserted=0)
    assert kept == []
    assert dropped[0].qualifiers["gene"] == "cds"


def test_an_empty_edit_does_not_drop_a_zero_length_feature() -> None:
    """`n^n+1` is a legitimate location to *load*; a no-op edit keeps it."""
    site = _ann(10, 10, "insertion_site")
    kept, dropped = shift_annotations([site], start=10, end=10, inserted=0)
    assert dropped == []
    assert (kept[0].start, kept[0].end) == (10, 10)


# ── a claim about the exact bases ───────────────────────────────────────────

def test_an_edit_inside_a_cds_does_not_keep_the_stale_translation() -> None:
    """Without the edited bases there is no way to restate it, so it goes."""
    kept, dropped = shift_annotations(
        [_cds()], start=SECOND_CODON[0], end=SECOND_CODON[0] + 1, inserted=0
    )
    assert kept == []
    assert dropped[0].qualifiers["translation"] == PROTEIN


def test_an_equal_length_replacement_inside_a_cds_is_not_silent() -> None:
    """delta == 0: nothing moves, and the record was still left wrong."""
    start, end = SECOND_CODON
    kept, dropped = shift_annotations([_cds()], start=start, end=end, inserted=end - start)
    assert kept == []
    assert dropped[0].feature_type == "CDS"


def test_given_the_edited_bases_the_translation_is_restated() -> None:
    start, end = SECOND_CODON
    edited = SEQ[:start] + "TGG" + SEQ[end:]
    kept, dropped = shift_annotations(
        [_cds()], start=start, end=end, inserted=3, new_sequence=edited
    )
    assert dropped == []
    assert kept[0].qualifiers["translation"] == "MW" + "A" * 9
    assert (kept[0].start, kept[0].end) == (CDS_START, CDS_END)


def test_a_restated_translation_is_what_the_bases_actually_encode() -> None:
    """The oracle is Biopython on the edited string, not the subject."""
    edited = SEQ[:20] + SEQ[23:]  # an in-frame codon deletion inside the CDS
    kept, _ = shift_annotations([_cds()], start=20, end=23, inserted=0, new_sequence=edited)
    coding = edited[CDS_START:CDS_END - 3]  # three bases shorter than it was
    assert kept[0].qualifiers["translation"] == str(
        Seq(coding).translate(table=1, to_stop=True)
    )
    assert kept[0].qualifiers["translation"] == "M" + "A" * 9


def test_a_minus_strand_spliced_cds_is_restated_in_reading_order() -> None:
    """Parts are stored 5'-to-3' (#93); the protein must follow that order."""
    rc = str(Seq(CDS_BASES).reverse_complement())
    seq = "CCCCC" + rc[:18] + "TTTT" + rc[18:] + "GGGGG"
    # Reading 3'-to-5' on the plus strand: the later interval comes first.
    parts = [(27, 27 + 18), (5, 23)]
    cds = _ann(5, 45, "orf", parts=parts, strand=Strand.MINUS,
               feature_type="CDS", translation=PROTEIN, codon_start="1")
    kept, dropped = shift_annotations(
        [cds], start=10, end=13, inserted=3, new_sequence=seq[:10] + "CCA" + seq[13:]
    )
    assert dropped == []
    edited = seq[:10] + "CCA" + seq[13:]
    bases = str(Seq(edited[27:45]).reverse_complement()) + str(
        Seq(edited[5:23]).reverse_complement()
    )
    assert kept[0].qualifiers["translation"] == str(
        Seq(bases).translate(table=1, to_stop=True)
    )


def test_a_non_standard_genetic_code_is_honoured() -> None:
    """/transl_table 2 reads TGA as tryptophan, not as a stop."""
    bases = "ATG" + "TGA" + "GCT" * 8 + "TAA"
    seq = "CCCCC" + bases + "GGGGG"
    cds = _ann(5, 5 + len(bases), "orf", feature_type="CDS",
               translation="MW" + "A" * 8, codon_start="1", transl_table="2")
    edited = seq[:11] + "GGT" + seq[14:]
    kept, dropped = shift_annotations([cds], start=11, end=14, inserted=3, new_sequence=edited)
    assert dropped == []
    assert kept[0].qualifiers["translation"] == "MWG" + "A" * 7


# ── what must not start happening ───────────────────────────────────────────

def test_an_insert_at_a_cds_boundary_leaves_it_alone() -> None:
    """Flush against a boundary the inserted bases are outside the feature."""
    for position in (CDS_START, CDS_END):
        kept, dropped = shift_annotations([_cds()], start=position, end=position, inserted=6)
        assert dropped == [], position
        assert kept[0].qualifiers["translation"] == PROTEIN, position


def test_a_feature_with_no_base_level_claim_still_grows_around_an_insert() -> None:
    """A repeat_region or a gene has no protein to invalidate: it survives."""
    kept, dropped = shift_annotations(
        [_ann(0, 20, "plain"), _ann(0, 20, "rep", feature_type="repeat_region")],
        start=10, end=10, inserted=5,
    )
    assert dropped == []
    assert [(a.start, a.end) for a in kept] == [(0, 25), (0, 25)]


def test_a_feature_elsewhere_still_shifts_by_exactly_delta() -> None:
    kept, dropped = shift_annotations(
        [_cds(), _ann(45, 50, "after")], start=0, end=3, inserted=0, new_sequence=SEQ[3:]
    )
    assert dropped == []
    assert [(a.start, a.end) for a in kept] == [(CDS_START - 3, CDS_END - 3), (42, 47)]
    assert kept[0].qualifiers["translation"] == PROTEIN
