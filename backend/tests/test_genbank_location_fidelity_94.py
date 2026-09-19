"""#94: a GenBank round-trip kept a location's integers and dropped the rest.

`_to_sequence` flattened every location to `(int, int)` pairs and `write_genbank`
rebuilt a plain `FeatureLocation`, so three things a location string carries
*besides* its coordinates were gone before the value was ever stored:

  <108..1007                     -> 108..1007      (the `<`/`>` position class)
  order(11..20,51..60)           -> join(11..20,51..60)   (the operator)
  join(201..210,J00194.1:100..202) -> join(201..210,100..202)  (the remote ref)

The last one is the worst behaved: it does not merely lose information, it hands
this record's own bases 100..202 to the feature.

This is the export half of the same reader/writer pair as #57 — #57 kept the
parts, #94 keeps everything else about them.

The oracle is Bio.SeqIO re-reading our own output, plus the literal FEATURES
lines, because `BeforePosition(5) == ExactPosition(5)` is `True`: any check that
compares coordinates is blind to a lost partial boundary.
"""

from __future__ import annotations

from io import StringIO

from Bio import SeqIO
from Bio.SeqFeature import AfterPosition, BeforePosition

from backend.formats.genbank import read_genbank, write_genbank
from backend.services.annotations import flip_annotations

BASES = ("atggcctgtg ggcatttggc caatttaggc catggacgtg gcatcacgtg gcatcacgtt\n"
         "       61 cggccatgga cgtaggcatc acgtggcatc acgtggccgg gaattcgatc ctgcaggtca\n"
         "      121 tggacgtggc atcacgtggc atcacgtggc cgggaattcg atcctgcagg tcatggacgt\n"
         "      181 ggcatcacgt ggcatcacgt ggccgggaat tcgatcctgc aggtcatgga cgtggcatca\n"
         "      241 cgtggcatca cgtggccggg aattcgatcc tgcaggtcat ggacgtggca tcacgtggca")

PARTIAL = f"""LOCUS       PARTIAL                  300 bp    DNA     linear   SYN 01-JAN-2024
DEFINITION  partial boundaries on both strands.
ACCESSION   PARTIAL
FEATURES             Location/Qualifiers
     source          1..300
     CDS             <108..200
                     /gene="runs_off_the_front"
     gene            complement(21..>60)
                     /gene="runs_off_the_back"
     exon            complement(<101..150)
                     /gene="partial_low_on_minus"
ORIGIN
        1 {BASES}
//
"""

OPERATORS = f"""LOCUS       OPERATOR                 300 bp    DNA     linear   SYN 01-JAN-2024
DEFINITION  an order() and a join() in one record.
ACCESSION   OPERATOR
FEATURES             Location/Qualifiers
     source          1..300
     misc_feature    order(11..20,51..60)
                     /label="ordered"
     misc_feature    join(101..110,141..150)
                     /label="joined"
ORIGIN
        1 {BASES}
//
"""

REMOTE = f"""LOCUS       REMOTE                   300 bp    DNA     linear   SYN 01-JAN-2024
DEFINITION  a join with one part in another entry.
ACCESSION   REMOTE
FEATURES             Location/Qualifiers
     source          1..300
     misc_feature    join(201..210,J00194.1:100..202)
                     /label="remote"
ORIGIN
        1 {BASES}
//
"""


def _round_trip(text: str):
    """(re-read export, export text) — never our output checked against itself."""
    exported = write_genbank(read_genbank(StringIO(text)))
    return SeqIO.read(StringIO(exported), "genbank"), exported


def _location_lines(genbank_text: str) -> list[str]:
    return [
        line[5:].strip().split(None, 1)[1]
        for line in genbank_text.splitlines()
        if line.startswith(" " * 5) and line[5:6] != " " and len(line[5:].split()) > 1
    ]


def test_a_partial_start_is_still_partial_after_export() -> None:
    """`<108..200` must not come back as `108..200`: INSDC's `<` says the feature
    starts before the first sequenced base, which an exact coordinate denies."""
    after, exported = _round_trip(PARTIAL)
    cds = after.features[1]
    assert isinstance(cds.location.start, BeforePosition), str(cds.location)
    assert "<108..200" in _location_lines(exported)


def test_a_partial_end_on_the_minus_strand_is_still_partial_after_export() -> None:
    after, exported = _round_trip(PARTIAL)
    gene = after.features[2]
    assert isinstance(gene.location.end, AfterPosition), str(gene.location)
    assert "complement(21..>60)" in _location_lines(exported)


def test_a_partial_lower_boundary_on_the_minus_strand_survives() -> None:
    after, exported = _round_trip(PARTIAL)
    exon = after.features[3]
    assert isinstance(exon.location.start, BeforePosition), str(exon.location)
    assert "complement(<101..150)" in _location_lines(exported)


def test_the_partial_boundaries_are_recorded_on_the_model() -> None:
    """The loss was at read time, so the flags have to be in the model, not
    reconstructed at export from something the writer cannot know."""
    anns = read_genbank(StringIO(PARTIAL))[0].annotations
    assert [(d.start_class, d.end_class) for d in anns[1].details_per_part()] == [("before", "exact")]
    assert [(d.start_class, d.end_class) for d in anns[2].details_per_part()] == [("exact", "after")]
    assert [(d.start_class, d.end_class) for d in anns[3].details_per_part()] == [("before", "exact")]


def test_an_exact_location_gains_no_angle_bracket() -> None:
    _, exported = _round_trip(OPERATORS)
    assert "<" not in exported.split("FEATURES")[1].split("ORIGIN")[0]
    assert ">" not in exported.split("FEATURES")[1].split("ORIGIN")[0]


def test_order_is_not_exported_as_join() -> None:
    """INSDC: join asserts the parts form one contiguous sequence, order asserts
    only that they occur in this order. Rewriting one as the other is a different
    claim about the molecule."""
    after, exported = _round_trip(OPERATORS)
    assert after.features[1].location.operator == "order", _location_lines(exported)
    assert "order(11..20,51..60)" in _location_lines(exported)


def test_a_real_join_still_exports_as_join() -> None:
    """The control: order() must not be fixed by breaking join()."""
    after, exported = _round_trip(OPERATORS)
    assert after.features[2].location.operator == "join"
    assert "join(101..110,141..150)" in _location_lines(exported)


def test_a_remote_reference_is_not_exported_as_local_bases() -> None:
    """`J00194.1:100..202` names bases in another entry; a bare `100..202` hands
    the feature this record's own bases instead, with nothing marking the swap."""
    _, exported = _round_trip(REMOTE)
    locations = _location_lines(exported)
    assert "join(201..210,100..202)" not in locations, locations
    assert "join(201..210,J00194.1:100..202)" in locations


def test_the_remote_accession_is_recorded_on_the_model() -> None:
    ann = read_genbank(StringIO(REMOTE))[0].annotations[1]
    details = ann.details_per_part()
    assert [d.ref for d in details] == [None, "J00194.1"]
    assert tuple(ann.parts[0]) == (200, 210), "the local part is this record's own"


def test_the_new_fields_default_to_the_old_behaviour() -> None:
    """Annotations from a reader with no notion of partiality (GFF, FASTA, the
    cloning planner) carry no details, and must still export as plain joins."""
    from backend.models.sequence import Annotation, MoleculeType, Sequence, Strand

    seq = Sequence(
        id="PLAIN", name="PLAIN", seq="acgt" * 25, molecule_type=MoleculeType.DNA,
        annotations=[Annotation(feature_type="gene", start=4, end=40,
                                parts=[(4, 10), (29, 40)], strand=Strand.PLUS)],
    )
    assert seq.annotations[0].part_details == []
    exported = write_genbank([seq])
    assert "join(5..10,30..40)" in _location_lines(exported)


def test_reverse_complement_moves_a_partial_boundary_to_the_other_end() -> None:
    """Reflecting `[s, e)` makes the old upper boundary the new lower one, so a
    `>` on the end has to come back as a `<` on the start — otherwise the flags
    go stale against the coordinates they describe."""
    seq = read_genbank(StringIO(PARTIAL))[0]
    gene = next(a for a in seq.annotations if a.qualifiers.get("gene") == "runs_off_the_back")
    flipped = flip_annotations([gene], len(seq.seq))[0]
    assert [(d.start_class, d.end_class) for d in flipped.details_per_part()] == [("before", "exact")]
