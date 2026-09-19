"""CO-1, CO-3, CO-4 — what a GenBank file's 1-based locations become in memory.

A GenBank file counts bases from 1 and includes both ends; the app's
``Annotation`` counts from 0 and excludes the end. Every number in here is the
pinned conversion of a location string in a real record, cross-checked against
``Bio.SeqIO`` run on the same bytes and against a re-slice of the raw sequence.

Subject: ``backend/formats/genbank.py:_to_sequence``.
"""

from __future__ import annotations

import pytest

from conftest import bio_extract, bio_parts, bio_strand, requirement, spliced

CO_1 = requirement("CO-1")
CO_3 = requirement("CO-3")
CO_4 = requirement("CO-4")
EXPECTED = CO_1["parameters"]["expected"]


def _ann(seq, i):
    return seq.annotations[i]


def test_co_1_feature_count_and_topology(jcv):
    assert len(jcv.annotations) == CO_1["parameters"]["feature_count"]
    assert jcv.length == CO_1["parameters"]["record_length"]
    assert jcv.is_circular is CO_1["parameters"]["is_circular"]


@pytest.mark.parametrize("row", EXPECTED, ids=[f"{r['i']}-{r['type']}" for r in EXPECTED])
def test_co_1_every_location_converts_to_zero_based_half_open(jcv, row):
    """`a..b` in the file must arrive as the interval (a-1, b)."""
    ann = _ann(jcv, row["i"])
    assert ann.feature_type == row["type"]
    assert [list(p) for p in ann.parts] == row["parts"], (
        f"{row['file_location']} should be {row['parts']}"
    )
    assert ann.strand.value == row["strand"]


@pytest.mark.parametrize("row", EXPECTED, ids=[f"{r['i']}-{r['type']}" for r in EXPECTED])
def test_co_1_load_agrees_with_biopython_on_the_same_bytes(jcv, jcv_bio, row):
    """The pinned table is not the only oracle: Biopython parses the file too."""
    ann, feat = _ann(jcv, row["i"]), jcv_bio.features[row["i"]]
    assert [tuple(p) for p in ann.parts] == bio_parts(feat)
    assert ann.strand.value == bio_strand(feat)


@pytest.mark.parametrize("row", EXPECTED, ids=[f"{r['i']}-{r['type']}" for r in EXPECTED])
def test_co_1_intervals_reslice_the_raw_sequence_to_the_right_bases(jcv, jcv_bio, row):
    """Re-derive the feature's bases by index; they must be the record's own."""
    ann = _ann(jcv, row["i"])
    by_index = spliced(jcv.seq, [tuple(p) for p in ann.parts], ann.strand.value)
    assert by_index == bio_extract(jcv_bio.features[row["i"]], jcv_bio)


@pytest.mark.parametrize("case", CO_3["parameters"]["cases"],
                         ids=lambda c: c["input"].rsplit("/", 1)[-1])
def test_co_3_origin_spanning_join_arrives_with_both_parts(case, jcv, bfdv_gb_text):
    """INSDC writes an origin-crossing feature as join(high..end,1..low)."""
    from backend.formats.genbank import read_genbank

    seq = jcv if "NC_001699" in case["input"] else read_genbank(bfdv_gb_text)[0]
    ann = seq.annotations[case["feature_index"]]
    assert ann.feature_type == case["type"]
    assert [list(p) for p in ann.parts] == case["parts"], case["file_location"]
    assert [e - s for s, e in ann.parts] == case["part_lengths"]
    assert sum(e - s for s, e in ann.parts) == case["total_length"]
    assert seq.is_circular is True
    # The outer bounds span the whole molecule; only `parts` carries the truth.
    assert (ann.start, ann.end) == (0, case["record_length"])


def test_co_4_complement_join_parts_are_in_reading_order(jcv, jcv_bio):
    """A minus-strand join lists its exons 5'-to-3', i.e. descending."""
    p = CO_4["parameters"]
    ann = jcv.annotations[p["feature_index"]]
    assert ann.strand.value == p["strand"]
    assert [list(q) for q in ann.parts] == p["in_memory_parts"]
    assert [e - s for s, e in ann.parts] == p["exon_lengths_in_reading_order"]
    assert bio_parts(jcv_bio.features[p["feature_index"]]) == [tuple(q) for q in ann.parts]


def test_co_4_complement_join_reads_the_protein_coding_strand(jcv, jcv_bio):
    """Concatenate-then-revcomp by index must give the CDS Biopython extracts."""
    p = CO_4["parameters"]
    ann = jcv.annotations[p["feature_index"]]
    by_index = spliced(jcv.seq, [tuple(q) for q in ann.parts], ann.strand.value)
    assert len(by_index) == p["extract_length"]
    assert by_index[:60] == p["extract_first_60"]
    assert by_index == bio_extract(jcv_bio.features[p["feature_index"]], jcv_bio)


@pytest.mark.parametrize("i", CO_4["parameters"]["simple_complement_features"])
def test_co_4_simple_complement_features_load_as_minus_strand(jcv, jcv_bio, i):
    assert jcv.annotations[i].strand.value == "-" == bio_strand(jcv_bio.features[i])
