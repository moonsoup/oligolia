"""#96: the HTTP edit endpoint stored a bare sequence.

`POST /sequences/{id}/edit` returns the edited bases *and* stores the result
under a new id. The stored record was built from five fields — id, name,
description, seq, molecule_type — so `annotations` and `is_circular` fell back
to their model defaults: an empty list and `False`. Edit an annotated circular
plasmid through the API and download the GenBank and every feature was gone and
the molecule was linear, with nothing said about it.

#59 fixed this for the GUI, which moves the annotations with the sequence via
`shift_annotations` / `flip_annotations`. These tests ask the same of the HTTP
path, over the same bytes the #93 tests use: a 120 bp circular record whose
features are a spliced CDS and an origin-spanning `rep_origin`, so "the
annotations came along" is not satisfiable by copying a list of coordinates.

The oracle is a re-slice of the edited string by index, never the subject's own
`extracted_bases`.
"""

from __future__ import annotations

from io import StringIO

from Bio.Seq import Seq
from fastapi.testclient import TestClient

from backend.formats.genbank import read_genbank

from ._spliced_circular_93 import ORIGIN_SPANNING, SPLICED_CDS, SPLICED_CIRCULAR

SEQ_ID = "edit96"


def _record() -> dict:
    """The pinned circular record, as the JSON body `POST /sequences/` takes."""
    seq = read_genbank(StringIO(SPLICED_CIRCULAR))[0]
    assert seq.is_circular and len(seq.annotations) == 3
    return seq.model_dump(mode="json") | {"id": SEQ_ID}


def _stored(client: TestClient, new_id: str) -> dict:
    r = client.get(f"/sequences/{new_id}")
    assert r.status_code == 200, r.text
    return r.json()


def _spliced(seq: str, parts: list[list[int]], strand: str) -> str:
    """What a feature denotes, built straight out of the string (INSDC join
    order; each part reverse-complemented on the minus strand)."""
    pieces = [seq[s:e] for s, e in parts]
    if strand == "-":
        return "".join(str(Seq(p).reverse_complement()) for p in pieces)
    return "".join(pieces)


def _edit(client: TestClient, **body) -> dict:
    client.post("/sequences/", json=_record())
    r = client.post(f"/sequences/{SEQ_ID}/edit", json=body)
    assert r.status_code == 200, r.text
    return r.json()


def _by_type(record: dict, feature_type: str) -> dict:
    matches = [a for a in record["annotations"] if a["feature_type"] == feature_type]
    assert len(matches) == 1, f"{feature_type}: {len(matches)} matches"
    return matches[0]


# ── the stored record carries the features ──────────────────────────────────

def test_insert_keeps_every_annotation(client: TestClient) -> None:
    _edit(client, operation="insert", position=30, insert_seq="TTTTT")
    stored = _stored(client, f"{SEQ_ID}_insert")
    assert len(stored["annotations"]) == 3, stored["annotations"]


def test_delete_keeps_the_untouched_annotations(client: TestClient) -> None:
    # [30, 40) sits between the CDS's two exons and touches nothing.
    _edit(client, operation="delete", position=30, end_position=40)
    stored = _stored(client, f"{SEQ_ID}_delete")
    assert len(stored["annotations"]) == 3, stored["annotations"]


def test_replace_keeps_the_untouched_annotations(client: TestClient) -> None:
    _edit(client, operation="replace", position=30, end_position=40, replacement="AAAAA")
    stored = _stored(client, f"{SEQ_ID}_replace")
    assert len(stored["annotations"]) == 3, stored["annotations"]


def test_reverse_complement_keeps_every_annotation(client: TestClient) -> None:
    _edit(client, operation="reverse_complement")
    stored = _stored(client, f"{SEQ_ID}_reverse_complement")
    assert len(stored["annotations"]) == 3, stored["annotations"]


# ── and they were shifted, not just copied ──────────────────────────────────

def test_an_inserts_downstream_exon_moves_and_the_upstream_one_does_not(
    client: TestClient,
) -> None:
    result = _edit(client, operation="insert", position=30, insert_seq="TTTTT")
    stored = _stored(client, f"{SEQ_ID}_insert")
    cds = _by_type(stored, "CDS")
    # CDS join(11..25,51..86) -> parts [(10, 25), (50, 86)]: the first exon ends
    # before the insert and must not move; the second starts after it and must
    # move by exactly the 5 bases inserted.
    assert [tuple(p) for p in cds["parts"]] == [(10, 25), (55, 91)]
    assert _spliced(result["result_seq"], cds["parts"], cds["strand"]) == SPLICED_CDS


def test_an_inserts_origin_spanning_feature_still_reads_its_own_bases(
    client: TestClient,
) -> None:
    result = _edit(client, operation="insert", position=30, insert_seq="TTTTT")
    stored = _stored(client, f"{SEQ_ID}_insert")
    origin = _by_type(stored, "rep_origin")
    # join(115..120,1..6) wraps the origin: the wrapped part is upstream of the
    # edit and stays, the other is downstream and shifts.
    assert [tuple(p) for p in origin["parts"]] == [(119, 125), (0, 6)]
    assert _spliced(result["result_seq"], origin["parts"], origin["strand"]) == ORIGIN_SPANNING


def test_a_deletes_downstream_feature_shifts_back_by_the_deleted_length(
    client: TestClient,
) -> None:
    result = _edit(client, operation="delete", position=30, end_position=40)
    stored = _stored(client, f"{SEQ_ID}_delete")
    cds = _by_type(stored, "CDS")
    assert [tuple(p) for p in cds["parts"]] == [(10, 25), (40, 76)]
    assert _spliced(result["result_seq"], cds["parts"], cds["strand"]) == SPLICED_CDS


def test_reverse_complement_flips_the_features_onto_the_other_strand(
    client: TestClient,
) -> None:
    result = _edit(client, operation="reverse_complement")
    stored = _stored(client, f"{SEQ_ID}_reverse_complement")
    cds = _by_type(stored, "CDS")
    assert cds["strand"] == "-"
    # [s, e) on a 120 bp record becomes [120 - e, 120 - s), and the part list
    # keeps its reading order rather than being re-sorted (#93) — so the
    # reflected second exon is still listed second.
    assert [tuple(p) for p in cds["parts"]] == [(95, 110), (34, 70)]
    assert _spliced(result["result_seq"], cds["parts"], cds["strand"]) == SPLICED_CDS


# ── topology ────────────────────────────────────────────────────────────────

def test_the_edited_plasmid_is_still_circular(client: TestClient) -> None:
    for body, new_id in (
        ({"operation": "insert", "position": 30, "insert_seq": "TTTTT"}, "insert"),
        ({"operation": "delete", "position": 30, "end_position": 40}, "delete"),
        ({"operation": "replace", "position": 30, "end_position": 40,
          "replacement": "AAAAA"}, "replace"),
        ({"operation": "reverse_complement"}, "reverse_complement"),
    ):
        _edit(client, **body)
        stored = _stored(client, f"{SEQ_ID}_{new_id}")
        assert stored["is_circular"] is True, f"{new_id} linearised the plasmid"


# ── the end of the chain: the file the user saves ───────────────────────────

def test_the_downloaded_genbank_of_an_edited_record_still_has_features(
    client: TestClient,
) -> None:
    _edit(client, operation="insert", position=30, insert_seq="TTTTT")
    stored = _stored(client, f"{SEQ_ID}_insert")
    r = client.post("/files/download/genbank", json=[stored])
    assert r.status_code == 200, r.text
    text = r.text
    assert "CDS" in text and "rep_origin" in text, text
    assert "circular" in text.splitlines()[0]
    # The stored record's name is "… (insert)" — a LOCUS name cannot hold a
    # space, and export used to raise ValueError instead of writing the file.
    reread = read_genbank(StringIO(text))[0]
    assert len(reread.annotations) == 3


def test_a_name_with_a_space_can_still_be_exported(client: TestClient) -> None:
    """The guard under the test above, without the edit endpoint in the way."""
    from backend.formats.genbank import write_genbank
    from backend.models.sequence import Sequence

    text = write_genbank([Sequence(id="spaced", name="my plasmid v2", seq="ACGTACGT")])
    assert text.startswith("LOCUS       my_plasmid_v2")


# ── nothing was taken away ──────────────────────────────────────────────────

def test_the_returned_bases_are_unchanged_by_all_of_this(client: TestClient) -> None:
    """The control: the endpoint's own result is what it always was."""
    original = read_genbank(StringIO(SPLICED_CIRCULAR))[0].seq
    result = _edit(client, operation="insert", position=30, insert_seq="TTTTT")
    assert result["result_seq"] == original[:30] + "TTTTT" + original[30:]
    assert result["diff_start"] == 30 and result["diff_end"] == 35


def test_translate_still_yields_a_protein_with_no_annotations(client: TestClient) -> None:
    """A translation changes the coordinate space from bases to residues, so it
    carries no features and is not a circular molecule."""
    client.post("/sequences/", json=_record())
    r = client.post(f"/sequences/{SEQ_ID}/edit", json={"operation": "translate"})
    assert r.status_code == 200, r.text
    stored = _stored(client, f"{SEQ_ID}_translate")
    assert stored["molecule_type"] == "PROTEIN"
    assert stored["annotations"] == []
    assert stored["is_circular"] is False


def test_transcribe_keeps_the_features_it_does_not_move(client: TestClient) -> None:
    """T→U over the whole molecule: same length, same coordinates."""
    client.post("/sequences/", json=_record())
    r = client.post(f"/sequences/{SEQ_ID}/edit", json={"operation": "transcribe"})
    assert r.status_code == 200, r.text
    stored = _stored(client, f"{SEQ_ID}_transcribe")
    assert stored["molecule_type"] == "RNA"
    assert len(stored["annotations"]) == 3
    assert [tuple(p) for p in _by_type(stored, "CDS")["parts"]] == [(10, 25), (50, 86)]
