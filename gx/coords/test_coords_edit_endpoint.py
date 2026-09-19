"""CO-14 — the edit endpoint's own copy of the record.

``POST /sequences/{id}/edit`` returns the edited bases *and* stores the result
under a new id (``backend/routers/sequences.py:136-145``). That stored record is
what ``GET /sequences/{id}`` hands back and what ``POST /files/download/genbank``
would write out, so whatever it drops is dropped from the file the user saves.

#59 was fixed in the GUI, where ``_commit_edit`` calls ``shift_annotations``.
These tests ask whether the HTTP path the same app exposes did the same.

Subject: ``backend/routers/sequences.py:edit_sequence``.
"""

from __future__ import annotations

import pytest

from conftest import requirement, spliced

CO_14 = requirement("CO-14")
P = CO_14["parameters"]
EDITS = P["edits"]


@pytest.fixture
def store(jcv):
    """A clean in-memory store holding the pinned JC polyomavirus record."""
    from backend.routers import sequences as router

    router._store.clear()
    jcv.id = P["stored_id"]
    router.add_sequence(jcv)
    assert len(router._store[P["stored_id"]].annotations) == P["source_annotation_count"]
    yield router
    router._store.clear()


def _edit(router, edit):
    from backend.models.sequence import SequenceEditRequest

    req = SequenceEditRequest(
        operation=edit["operation"],
        position=edit.get("position"),
        end_position=edit.get("end_position"),
        insert_seq=edit.get("insert_seq"),
        replacement=edit.get("replacement"),
    )
    return router.edit_sequence(P["stored_id"], req)


@pytest.mark.parametrize("edit", EDITS, ids=lambda e: e["operation"])
def test_co_14_the_edit_itself_is_right(store, edit):
    """Control: the bases the endpoint returns are correct, so the failures
    below are about the annotations, not about the splice."""
    result = _edit(store, edit)
    assert len(result.result_seq) == edit["expected_length"]
    assert store._store[edit["new_id"]].seq == result.result_seq


@pytest.mark.parametrize("edit", EDITS, ids=lambda e: e["operation"])
def test_co_14_the_stored_record_keeps_its_annotations(store, edit):
    _edit(store, edit)
    stored = store._store[edit["new_id"]]
    assert len(stored.annotations) >= P["expected_stored_annotation_count_at_least"], (
        f"{edit['operation']} stored {edit['new_id']} with "
        f"{len(stored.annotations)} of {P['source_annotation_count']} features"
    )


@pytest.mark.parametrize("edit", EDITS, ids=lambda e: e["operation"])
def test_co_14_the_stored_record_keeps_its_topology(store, edit):
    _edit(store, edit)
    assert store._store[edit["new_id"]].is_circular is P["expected_is_circular"], (
        f"{edit['operation']} turned a circular plasmid into a linear one"
    )


def test_co_14_surviving_annotations_point_at_the_right_bases(store, jcv_bio):
    """If features are carried across an insert, they must be shifted: every
    feature downstream of the insert still has to read its own bases."""
    edit = EDITS[0]
    before = {
        (f.type, str(f.location)): str(f.extract(jcv_bio.seq)) for f in jcv_bio.features
    }
    _edit(store, edit)
    stored = store._store[edit["new_id"]]
    if not stored.annotations:
        pytest.fail("no annotations survived the edit; nothing to check")
    wrong = [
        a.feature_type
        for a in stored.annotations
        if spliced(stored.seq, [tuple(p) for p in a.parts], a.strand.value)
        not in before.values()
    ]
    assert wrong == [], wrong


def test_co_14_a_downloaded_genbank_of_the_edited_record_still_has_features(store):
    """The end of the chain: what the user actually saves."""
    from backend.formats.genbank import write_genbank

    from conftest import feature_lines

    edit = EDITS[0]
    _edit(store, edit)
    exported = write_genbank([store._store[edit["new_id"]]])
    assert feature_lines(exported), (
        "the exported GenBank file has an empty FEATURES block"
    )
