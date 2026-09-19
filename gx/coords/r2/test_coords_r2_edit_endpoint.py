"""CO2-1, CO2-2, CO2-3 — what ``POST /sequences/{id}/edit`` owes each feature.

Round 2. The sealed round-1 test ``gx-co-53`` asked every surviving annotation
to extract bases the *pre-edit* record already held, which no feature that
contains the edit can do and still be the feature it was. These three
requirements say instead what CO-14 actually says: a containing feature grows
with the edit, a downstream one shifts, an upstream one does not move — the
arithmetic ``backend/services/annotations.py:shift_annotations`` performs and
``gui/panels/sequence_panel.py:_commit_edit`` has called since #59.

Subject: ``backend/routers/sequences.py:edit_sequence``.
"""

from __future__ import annotations

import pytest
from coords_r2_oracle import (
    bio_extract,
    bio_parts,
    bio_strand,
    contains,
    expected_extract,
    expected_parts,
    jcv_bio,
    jcv_text,
    requirement,
    sha256,
    spliced,
)

CO2_1 = requirement("CO2-1")
CO2_2 = requirement("CO2-2")
CO2_3 = requirement("CO2-3")

INSERT = CO2_1["parameters"]["edit"]
CONTAINING = CO2_1["parameters"]["containing_features"]
EDITS = CO2_2["parameters"]["edits"]
PER_EDIT = {e["name"]: e for e in CO2_2["parameters"]["per_edit"]}

BIO = jcv_bio()
SEQ = str(BIO.seq)


# ── fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def store():
    """A clean in-memory store holding the app's parse of the pinned record."""
    from backend.formats.genbank import read_genbank
    from backend.routers import sequences as router

    seqs = read_genbank(jcv_text())
    assert len(seqs) == 1
    seq = seqs[0]
    seq.id = CO2_1["parameters"]["stored_id"]
    router._store.clear()
    router.add_sequence(seq)
    assert (
        len(router._store[seq.id].annotations)
        == CO2_1["parameters"]["source_annotation_count"]
    )
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
    router.edit_sequence(CO2_1["parameters"]["stored_id"], req)
    return router._store[edit["new_id"]]


def _stored_in_file_order(router, edit):
    """The edited record, with its annotations indexable by original position.

    None of the three pinned edits partially overlaps any feature of
    NC_001699.1, so no feature has grounds to be dropped and stored index *i* is
    original index *i*. The assertion makes that mapping a checked precondition
    rather than an assumption: if it ever fails, every index below is read from
    the register, not guessed.
    """
    stored = _edit(router, edit)
    assert len(stored.annotations) == len(BIO.features), (
        f"{edit['name']} stored {len(stored.annotations)} of {len(BIO.features)} "
        "features; no feature of this record partially overlaps any pinned edit, "
        "so nothing had grounds to be dropped"
    )
    return stored


def _stored_parts(ann):
    return [tuple(p) for p in ann.parts]


# ── CO2-1 — a feature containing the insert grows by the inserted bases ──────

def test_co2_1_the_insert_keeps_every_feature_and_the_topology(store):
    """Nothing is dropped: the five features that contain the insert are kept,
    growing, not discarded. Dropping them would lose the whole-plasmid `source`
    feature on every insertion — the loss #96 exists to fix."""
    p = CO2_1["parameters"]
    stored = _edit(store, INSERT)
    assert len(stored.seq) == INSERT["expected_length"]
    assert len(stored.annotations) == p["expected_stored_annotation_count"]
    assert stored.is_circular is p["expected_is_circular"]
    kept_types = [a.feature_type for a in stored.annotations]
    for row in CONTAINING:
        assert kept_types[row["index"]] == row["feature_type"], (
            f"feature {row['index']} ({row['feature_type']}) — which strictly "
            f"contains position {INSERT['position']} — is not at its own index"
        )


@pytest.mark.parametrize(
    "row", CONTAINING, ids=lambda r: f"{r['index']}-{r['feature_type']}"
)
def test_co2_1_a_containing_feature_grows_by_the_inserted_bases(store, row):
    """The containing part keeps its start and moves its end by +30; any other
    part of the same feature follows the ordinary downstream/upstream rule."""
    before = bio_parts(BIO.features[row["index"]])
    assert [list(p) for p in before] == row["parts_before"]  # register vs Biopython
    want = expected_parts(before, INSERT["start"], INSERT["end"], INSERT["delta"])
    assert [list(p) for p in want] == row["parts_after"]  # register vs re-derivation

    stored = _stored_in_file_order(store, INSERT)
    got = _stored_parts(stored.annotations[row["index"]])
    assert got == want, (
        f"{row['feature_type']} {before} contains position {INSERT['position']}; "
        f"expected {want}, got {got}"
    )
    assert any(contains(p, INSERT["start"], INSERT["end"]) for p in before)


@pytest.mark.parametrize(
    "row", CONTAINING, ids=lambda r: f"{r['index']}-{r['feature_type']}"
)
def test_co2_1_a_containing_feature_still_describes_its_own_feature(store, row):
    """Its bases must be its own pre-edit bases with the insert spliced in at the
    corresponding offset — 30 longer, and byte-identical nowhere else."""
    feature = BIO.features[row["index"]]
    before = bio_extract(feature, BIO)
    want = expected_extract(
        bio_parts(feature), bio_strand(feature), SEQ,
        INSERT["start"], INSERT["end"], INSERT["inserted_bases"],
    )
    # The register's pinned digests and the re-derivation must agree before the
    # subject is consulted at all.
    assert sha256(before) == row["extract_sha256_before"]
    assert sha256(want) == row["extract_sha256_after"]
    assert (len(before), len(want)) == (
        row["extract_length_before"], row["extract_length_after"],
    )
    assert len(want) == len(before) + INSERT["delta"]

    stored = _stored_in_file_order(store, INSERT)
    ann = stored.annotations[row["index"]]
    got = spliced(stored.seq, _stored_parts(ann), ann.strand.value)
    assert got == want, (
        f"{row['feature_type']} reads {len(got)} bases after the insert, "
        f"expected {len(want)} (its own {len(before)} with the 30 inserted bases "
        "spliced in)"
    )


# ── CO2-2 — a part strictly downstream shifts by delta ───────────────────────

@pytest.mark.parametrize("edit", EDITS, ids=lambda e: e["name"])
def test_co2_2_strictly_downstream_parts_shift_by_exactly_delta(store, edit):
    """Re-derive each shifted interval by index: new[s+d:e+d] == old[s:e]."""
    stored = _stored_in_file_order(store, edit)
    delta, end = edit["delta"], edit["end"]
    wrong = []
    for i, feature in enumerate(BIO.features):
        got = _stored_parts(stored.annotations[i])
        for (p_start, p_end), (g_start, g_end) in zip(bio_parts(feature), got):
            if p_start < end:
                continue  # not strictly downstream; CO2-1/CO2-3's business
            if (g_start, g_end) != (p_start + delta, p_end + delta):
                wrong.append((i, feature.type, (p_start, p_end), (g_start, g_end)))
            elif stored.seq[g_start:g_end] != SEQ[p_start:p_end]:
                wrong.append((i, feature.type, "bases differ at the shifted interval"))
    assert wrong == [], wrong[:5]


@pytest.mark.parametrize("edit", EDITS, ids=lambda e: e["name"])
def test_co2_2_a_feature_that_does_not_contain_the_edit_extracts_byte_identically(
    store, edit
):
    """No part of it changed, so nothing about the bases it reads may change."""
    stored = _stored_in_file_order(store, edit)
    per = PER_EDIT[edit["name"]]
    untouched = set(per["strictly_downstream_feature_indices"]) | set(
        per["strictly_upstream_feature_indices"]
    ) | set(per["straddling_feature_indices"])
    assert untouched.isdisjoint(per["containing_feature_indices"])

    wrong = []
    for i in sorted(untouched):
        ann = stored.annotations[i]
        got = spliced(stored.seq, _stored_parts(ann), ann.strand.value)
        want = bio_extract(BIO.features[i], BIO)
        if got != want:
            wrong.append((i, BIO.features[i].type, len(want), len(got)))
    assert wrong == [], f"{len(wrong)} features changed bases: {wrong[:5]}"


def test_co2_2_worked_example_the_origin_spanning_join_moves_only_one_part(store):
    """``join(5118..5130,1..12)``: the first part is downstream of the insert and
    moves by +30, the second is upstream and does not — and the feature still
    reads exactly the same 25 bases."""
    ex = CO2_2["parameters"]["worked_example_origin_spanning_join"]
    feature = BIO.features[ex["index"]]
    assert feature.type == ex["feature_type"]
    assert [list(p) for p in bio_parts(feature)] == ex["parts_before"]

    stored = _stored_in_file_order(store, INSERT)
    ann = stored.annotations[ex["index"]]
    assert [list(p) for p in _stored_parts(ann)] == ex["parts_after_insert_30_at_3000"]
    got = spliced(stored.seq, _stored_parts(ann), ann.strand.value)
    assert sha256(got) == ex["extract_sha256"]
    assert got == bio_extract(feature, BIO)


# ── CO2-3 — a feature strictly upstream is untouched ─────────────────────────

@pytest.mark.parametrize("edit", EDITS, ids=lambda e: e["name"])
def test_co2_3_strictly_upstream_features_are_untouched(store, edit):
    """Same parts, same strand, same bases. An edit downstream of a feature is
    none of that feature's business."""
    stored = _stored_in_file_order(store, edit)
    indices = PER_EDIT[edit["name"]]["strictly_upstream_feature_indices"]
    assert indices, "the register pins no upstream feature for this edit"

    wrong = []
    for i in indices:
        feature = BIO.features[i]
        want_parts = bio_parts(feature)
        assert all(p_end <= edit["start"] for _s, p_end in want_parts)  # fixture sanity
        ann = stored.annotations[i]
        got_parts = _stored_parts(ann)
        got = spliced(stored.seq, got_parts, ann.strand.value)
        if (got_parts, ann.strand.value, got) != (
            want_parts, bio_strand(feature), bio_extract(feature, BIO)
        ):
            wrong.append((i, feature.type, want_parts, got_parts))
    assert wrong == [], f"{len(wrong)} upstream features moved: {wrong[:5]}"
