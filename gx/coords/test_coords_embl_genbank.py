"""CO-13 — one INSDC record, two flat-file formats, one molecule.

AF071878.1 as NCBI serves it (GenBank) and as EBI serves it (EMBL) are two
independent renderings of the same submission. The subject has a separate
reader for each; they must agree about the sequence, the topology and every
feature the two archives share.

The two archives do differ: NCBI's copy carries `gene` features ENA's does not.
That is a real content difference between databases, so the comparison is over
the feature keys both files contain — which ``Bio.SeqIO`` identifies on each
file independently, before the subject is involved.

Subject: ``backend/formats/genbank.py:read_genbank`` / ``read_embl`` (#67).
"""

from __future__ import annotations

import pytest

from conftest import (
    BFDV_EMBL,
    BFDV_GB,
    bio_extract,
    bio_parts,
    bio_record,
    bio_strand,
    requirement,
    spliced,
)

CO_13 = requirement("CO-13")
P = CO_13["parameters"]
SHARED = P["shared_features"]


@pytest.fixture(scope="module")
def pair():
    from backend.formats.genbank import read_embl, read_genbank

    return read_genbank(BFDV_GB.read_text())[0], read_embl(BFDV_EMBL.read_text())[0]


def test_co_13_the_two_files_really_are_the_same_record():
    """Fixture sanity, from the authority: same bases, same topology."""
    gb, em = bio_record(BFDV_GB), bio_record(BFDV_EMBL, "embl")
    assert str(gb.seq) == str(em.seq)
    assert len(gb.seq) == P["record_length"]
    assert gb.annotations.get("topology") == em.annotations.get("topology") == "circular"
    assert (len(gb.features), len(em.features)) == (
        P["genbank_feature_count"], P["embl_feature_count"]
    )


def test_co_13_both_readers_return_one_record_of_the_same_sequence(pair):
    gb, em = pair
    assert gb.seq == em.seq
    assert gb.length == em.length == P["record_length"]


def test_co_13_both_readers_agree_the_record_is_circular(pair):
    gb, em = pair
    assert gb.is_circular is em.is_circular is P["is_circular"]


@pytest.mark.parametrize("row", SHARED, ids=lambda r: f"{r['type']}-{r['file_location']}")
def test_co_13_shared_features_load_identically_from_both_formats(pair, row):
    gb, em = pair
    from_gb = [a for a in gb.annotations
               if a.feature_type == row["type"] and [list(p) for p in a.parts] == row["parts"]]
    from_em = [a for a in em.annotations
               if a.feature_type == row["type"] and [list(p) for p in a.parts] == row["parts"]]
    assert from_gb, f"GenBank reader lost {row['type']} {row['file_location']}"
    assert from_em, f"EMBL reader lost {row['type']} {row['file_location']}"
    assert from_gb[0].strand.value == row["strand"]
    assert from_em[0].strand.value == row["strand"]


@pytest.mark.parametrize("row", SHARED, ids=lambda r: f"{r['type']}-{r['file_location']}")
def test_co_13_shared_features_denote_the_same_bases_in_both_formats(pair, row):
    gb, em = pair
    parts = [tuple(p) for p in row["parts"]]
    from_gb = spliced(gb.seq, parts, row["strand"])
    from_em = spliced(em.seq, parts, row["strand"])
    assert from_gb == from_em

    bio = bio_record(BFDV_GB)
    feat = next(f for f in bio.features
                if f.type == row["type"] and bio_parts(f) == parts
                and bio_strand(f) == row["strand"])
    assert from_gb == bio_extract(feat, bio)


def test_co_13_the_embl_reader_matches_biopython_on_the_embl_bytes(pair):
    """Independently of the GenBank copy: the .embl file's own oracle."""
    _gb, em = pair
    bio = bio_record(BFDV_EMBL, "embl")
    assert [(a.feature_type, [tuple(p) for p in a.parts], a.strand.value)
            for a in em.annotations] == [
        (f.type, bio_parts(f), bio_strand(f)) for f in bio.features
    ]


def test_co_13_the_origin_spanning_join_survives_the_embl_reader_too(pair):
    _gb, em = pair
    row = next(r for r in SHARED if r["type"] == "stem_loop")
    ann = next(a for a in em.annotations if a.feature_type == "stem_loop")
    assert [list(p) for p in ann.parts] == row["parts"]
    assert em.is_circular is True
