"""#57: a GenBank round-trip silently corrupted feature locations.

`_to_sequence` read only `location.start`/`location.end`, so every `join()`
collapsed to its outer bounds:

  join(5..10,30..40)     -> 5..40, swallowing the intron
  join(91..100,1..20)    -> 1..100 on a 100 bp circular record, i.e. the whole
                            plasmid, which the plasmid map then drew as a full circle

and the writer stringified multi-value qualifiers, so two /note entries came back
as one /note="['first', 'second']".

The oracle throughout is Bio.SeqIO: our output is re-read with Biopython and the
parts must survive.
"""

from __future__ import annotations

from io import StringIO

from Bio import SeqIO

from backend.formats.genbank import read_genbank, write_genbank

SPLICED = """LOCUS       SPLICED                   60 bp    DNA     linear   SYN 01-JAN-2024
DEFINITION  a feature with an intron.
ACCESSION   SPLICED
FEATURES             Location/Qualifiers
     source          1..60
     CDS             join(5..10,30..40)
                     /gene="demo"
                     /note="first"
                     /note="second"
ORIGIN
        1 atggcctgtg ggcatttggc caatttaggc catggacgtg gcatcacgtg gcatcacgtt
//
"""

WRAPPED = """LOCUS       PLASMID                  100 bp    DNA     circular SYN 01-JAN-2024
DEFINITION  a feature spanning the origin.
ACCESSION   PLASMID
FEATURES             Location/Qualifiers
     source          1..100
     gene            join(91..100,1..20)
                     /gene="wraps"
ORIGIN
        1 atggcctgtg ggcatttggc caatttaggc catggacgtg gcatcacgtg gcatcacgtt
       61 cggccatgga cgtaggcatc acgtggcatc acgtggccgg
//
"""


def _feature(seq, ftype: str):
    return next(a for a in seq.annotations if a.feature_type == ftype)


def test_a_join_keeps_both_parts(client=None) -> None:
    seqs = read_genbank(StringIO(SPLICED))
    cds = _feature(seqs[0], "CDS")

    # 0-based half-open, as Biopython reports them.
    assert cds.parts == [(4, 10), (29, 40)], cds.parts


def test_the_outer_bounds_are_unchanged_for_existing_consumers() -> None:
    """start/end must keep meaning what they did, so the plasmid map and the
    feature table are unaffected by the new field."""
    cds = _feature(read_genbank(StringIO(SPLICED))[0], "CDS")
    assert cds.start == 4
    assert cds.end == 40


def test_a_simple_feature_reports_one_part() -> None:
    src = _feature(read_genbank(StringIO(SPLICED))[0], "source")
    assert src.parts == [(0, 60)], src.parts


def test_a_join_survives_export_and_reimport() -> None:
    """The invariant: load -> export -> reload must not change the parts."""
    original = _feature(read_genbank(StringIO(SPLICED))[0], "CDS")
    text = write_genbank(read_genbank(StringIO(SPLICED)))
    reloaded = _feature(read_genbank(StringIO(text))[0], "CDS")
    assert reloaded.parts == original.parts, (reloaded.parts, original.parts)


def test_biopython_sees_the_join_in_our_output() -> None:
    """Differential check: the exported text must really contain a join()."""
    text = write_genbank(read_genbank(StringIO(SPLICED)))
    assert "join(" in text, text

    rec = next(SeqIO.parse(StringIO(text), "genbank"))
    cds = next(f for f in rec.features if f.type == "CDS")
    parts = [(int(p.start), int(p.end)) for p in cds.location.parts]
    assert parts == [(4, 10), (29, 40)], parts


def test_an_origin_spanning_feature_is_not_the_whole_plasmid() -> None:
    """#57's worst case: join(91..100,1..20) loaded as [0,100) — everything."""
    seq = read_genbank(StringIO(WRAPPED))[0]
    gene = _feature(seq, "gene")

    assert gene.parts == [(90, 100), (0, 20)], gene.parts
    covered = sum(e - s for s, e in gene.parts)
    assert covered == 30, f"the feature covers 30 bp, not {covered}"
    assert covered < len(seq.seq), "a wrap-around feature is not the entire record"


def test_an_origin_spanning_feature_survives_the_round_trip() -> None:
    original = _feature(read_genbank(StringIO(WRAPPED))[0], "gene")
    text = write_genbank(read_genbank(StringIO(WRAPPED)))
    reloaded = _feature(read_genbank(StringIO(text))[0], "gene")
    assert reloaded.parts == original.parts, (reloaded.parts, original.parts)


def test_multi_value_qualifiers_stay_multi_value() -> None:
    """#57: two /note entries re-exported as /note="['first', 'second']"."""
    cds = _feature(read_genbank(StringIO(SPLICED))[0], "CDS")
    assert cds.qualifiers["note"] == ["first", "second"], cds.qualifiers["note"]

    text = write_genbank(read_genbank(StringIO(SPLICED)))
    assert "['first'" not in text, "a qualifier list was stringified into one value"
    assert text.count('/note="first"') == 1, text
    assert text.count('/note="second"') == 1, text

    reloaded = _feature(read_genbank(StringIO(text))[0], "CDS")
    assert reloaded.qualifiers["note"] == ["first", "second"], reloaded.qualifiers


def test_a_single_value_qualifier_stays_scalar() -> None:
    """Existing behaviour: a lone value is unwrapped for convenience."""
    cds = _feature(read_genbank(StringIO(SPLICED))[0], "CDS")
    assert cds.qualifiers["gene"] == "demo", cds.qualifiers["gene"]
