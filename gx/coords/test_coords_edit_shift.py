"""CO-8, CO-9, CO-10, CO-11 — what an edit does to the annotations around it.

#59 fixed the case where nothing moved at all. These tests ask the next
question: after the move, does every feature still describe the bases it
claims? The oracle is the *unedited* record's extracted subsequences (from
``Bio.SeqIO`` on the pinned bytes) and a re-slice of the edited string by
index — the subject's own reported coordinates are never the standard.

Subject: ``backend/services/annotations.py:shift_annotations``, which
``gui/panels/sequence_panel.py:_commit_edit`` calls on every insert, delete and
replace.
"""

from __future__ import annotations

import pytest
from Bio.Seq import Seq

from conftest import (
    bio_extract,
    edited_sequence,
    feature_lines,
    requirement,
    spliced,
    touched,
)

CO_8 = requirement("CO-8")
CO_9 = requirement("CO-9")
CO_10 = requirement("CO-10")
CO_11 = requirement("CO-11")

EDITS = CO_8["parameters"]["edits"]
EDIT_IDS = [e["name"] for e in EDITS]


def _apply(seq, edit):
    """Run one pinned edit, and report which original index each kept feature is."""
    new, dropped = edited_sequence(
        seq,
        start=edit["start"],
        end=edit["end"],
        inserted=edit.get("inserted_bases", ""),
    )
    dropped_ids = {id(a) for a in dropped}
    origin = [i for i, a in enumerate(seq.annotations) if id(a) not in dropped_ids]
    return new, dropped, origin


@pytest.mark.parametrize("edit", EDITS, ids=EDIT_IDS)
def test_co_8_untouched_features_keep_byte_identical_bases(jcv, jcv_bio, edit):
    """An edit elsewhere in the molecule must not change what a feature reads."""
    before = [bio_extract(f, jcv_bio) for f in jcv_bio.features]
    parts_before = [[tuple(p) for p in a.parts] for a in jcv.annotations]
    new, dropped, origin = _apply(jcv, edit)

    wrong = []
    for j, ann in enumerate(new.annotations):
        i = origin[j]
        if touched(parts_before[i], edit["start"], edit["end"]):
            continue
        after = spliced(new.seq, [tuple(p) for p in ann.parts], ann.strand.value)
        if after != before[i]:
            wrong.append((i, ann.feature_type, len(before[i]), len(after)))
    assert wrong == [], f"{len(wrong)} untouched features changed bases: {wrong[:5]}"


@pytest.mark.parametrize("edit", EDITS, ids=EDIT_IDS)
def test_co_8_an_untouched_feature_is_never_dropped(jcv, edit):
    parts_before = [[tuple(p) for p in a.parts] for a in jcv.annotations]
    _, dropped, _ = _apply(jcv, edit)
    ids = {id(a) for a in dropped}
    surprises = [
        (i, a.feature_type)
        for i, a in enumerate(jcv.annotations)
        if id(a) in ids and not touched(parts_before[i], edit["start"], edit["end"])
    ]
    assert surprises == []


@pytest.mark.parametrize("edit", EDITS, ids=EDIT_IDS)
def test_co_9_downstream_parts_shift_by_exactly_delta(jcv, edit):
    """Re-derive each shifted interval by index: new[s+d:e+d] == old[s:e]."""
    delta, start, end = edit["delta"], edit["start"], edit["end"]
    parts_before = [[tuple(p) for p in a.parts] for a in jcv.annotations]
    old_seq = jcv.seq
    new, _, origin = _apply(jcv, edit)

    wrong = []
    for j, ann in enumerate(new.annotations):
        i = origin[j]
        got = [tuple(p) for p in ann.parts]
        for (p_start, p_end), (g_start, g_end) in zip(parts_before[i], got):
            if p_start >= end:
                want = (p_start + delta, p_end + delta)
            elif p_end <= start:
                want = (p_start, p_end)
            else:
                continue  # a part containing the edit is CO-11's business
            if (g_start, g_end) != want:
                wrong.append((i, ann.feature_type, (p_start, p_end), (g_start, g_end), want))
            elif new.seq[g_start:g_end] != old_seq[p_start:p_end]:
                wrong.append((i, ann.feature_type, "bases differ at the shifted interval"))
    assert wrong == [], wrong[:5]


def test_co_9_worked_example_insert_upstream_of_an_origin_spanning_join(jcv):
    """Only the part downstream of the insert moves; the wrapped part does not."""
    ex = CO_9["parameters"]["worked_example"]
    edit = next(e for e in EDITS if e["name"] == ex["edit"])
    new, _, _ = _apply(jcv, edit)
    ann = new.annotations[ex["feature_index"]]
    assert [list(p) for p in ann.parts] == ex["after"]


def test_co_9_worked_example_delete_upstream_of_a_simple_feature(jcv):
    from backend.formats.genbank import write_genbank

    ex = CO_9["parameters"]["second_worked_example"]
    edit = next(e for e in EDITS if e["name"] == ex["edit"])
    new, _, origin = _apply(jcv, edit)
    j = origin.index(ex["feature_index"])
    assert [list(p) for p in new.annotations[j].parts] == ex["after"]
    assert feature_lines(write_genbank([new]))[j][1] == ex["exported_after"]


@pytest.mark.parametrize(
    "case",
    [
        {"name": "exact_delete_of_repeat_region_5074..5090",
         "start": CO_10["parameters"]["edit"]["start"],
         "end": CO_10["parameters"]["edit"]["end"],
         "feature_index": CO_10["parameters"]["feature_index"],
         "forbidden_exported_location": CO_10["parameters"]["forbidden_exported_location"]},
        *CO_10["parameters"]["companion_cases"],
    ],
    ids=lambda c: c["name"],
)
def test_co_10_a_wholly_deleted_feature_is_reported_not_kept_empty(jcv, case):
    """Delete exactly a feature's bases and the feature has nothing left to
    describe; it must come back in the dropped list."""
    edit = {"start": case["start"], "end": case["end"], "inserted_bases": ""}
    _, dropped, _ = _apply(jcv, edit)
    target = jcv.annotations[case["feature_index"]]
    assert any(a is target for a in dropped), (
        f"{target.feature_type} {[tuple(p) for p in target.parts]} survived a "
        f"delete of exactly its own bases [{case['start']}, {case['end']})"
    )


@pytest.mark.parametrize(
    "case",
    [
        {"name": "exact_delete_of_repeat_region_5074..5090",
         "start": CO_10["parameters"]["edit"]["start"],
         "end": CO_10["parameters"]["edit"]["end"],
         "forbidden_exported_location": CO_10["parameters"]["forbidden_exported_location"]},
        *CO_10["parameters"]["companion_cases"],
    ],
    ids=lambda c: c["name"],
)
def test_co_10_no_zero_length_between_position_is_exported(jcv, case):
    """INSDC's `n^n+1` points at a site *between* two bases. A feature that was
    a range must never be written that way."""
    from backend.formats.genbank import write_genbank

    edit = {"start": case["start"], "end": case["end"], "inserted_bases": ""}
    new, _, _ = _apply(jcv, edit)
    locs = [loc for _key, loc in feature_lines(write_genbank([new]))]
    assert case["forbidden_exported_location"] not in locs, locs


def test_co_10_no_exported_location_names_base_zero(jcv):
    """A GenBank file counts from 1: `0^1` is not a location in any record."""
    from backend.formats.genbank import write_genbank

    case = CO_10["parameters"]["companion_cases"][1]  # delete of rep_origin 1..12
    new, _, _ = _apply(jcv, {"start": case["start"], "end": case["end"], "inserted_bases": ""})
    offending = [loc for _k, loc in feature_lines(write_genbank([new]))
                 if "0^1" in loc or loc.startswith("0.")]
    assert offending == [], offending


def test_co_11_a_cds_whose_bases_changed_does_not_keep_a_stale_translation(jcv):
    """Either drop the CDS or restate its protein — do not export a
    /translation that no longer translates from the bases it sits on."""
    from backend.formats.genbank import write_genbank

    from conftest import reread

    p = CO_11["parameters"]
    edit = {"start": p["edit"]["start"], "end": p["edit"]["end"], "inserted_bases": ""}
    target = jcv.annotations[p["feature_index"]]
    assert target.feature_type == p["feature_type"]
    assert str(target.qualifiers["translation"]).startswith(p["declared_translation_prefix"])

    new, dropped, origin = _apply(jcv, edit)
    if any(a is target for a in dropped):
        return  # reported: acceptable

    j = origin.index(p["feature_index"])
    exported = reread(write_genbank([new]))
    feat = exported.features[j]
    declared = feat.qualifiers["translation"][0]
    actual = str(feat.extract(exported.seq).translate(table=p["translation_table"], to_stop=True))
    assert declared == actual, (
        f"CDS kept as {feat.location} ({len(feat.extract(exported.seq))} bases) with "
        f"the pre-edit /translation; bases now encode {actual[:24]}…"
    )


def test_co_11_the_shortened_cds_is_not_left_out_of_frame(jcv):
    """A one-base deletion inside a CDS leaves a length no codon reading can
    use; keeping it silently is the defect, whichever way it is fixed."""
    p = CO_11["parameters"]
    target = jcv.annotations[p["feature_index"]]
    # Only fair for a complete CDS: a partial one (`<`/`>`) may legitimately
    # not be a whole number of codons. This one is neither partial nor short.
    assert "<" not in str(target.qualifiers) and p["length_before"] % 3 == 0
    edit = {"start": p["edit"]["start"], "end": p["edit"]["end"], "inserted_bases": ""}
    new, dropped, origin = _apply(jcv, edit)
    if any(a is target for a in dropped):
        return
    ann = new.annotations[origin.index(p["feature_index"])]
    length = sum(e - s for s, e in ann.parts)
    assert length == p["length_after_if_adjusted"]  # fixture sanity
    assert length % 3 == 0, (
        f"CDS kept at {length} bases, not a multiple of three, with its original "
        "/codon_start and /translation"
    )


def test_co_11_a_same_length_replace_inside_a_cds_is_not_silent(jcv):
    """delta == 0, so nothing moves — but the CDS's own bases changed."""
    from backend.formats.genbank import write_genbank

    from conftest import reread

    c = CO_11["parameters"]["companion_case"]
    replacement = "A" * (c["end"] - c["start"])
    new, dropped, origin = _apply(
        jcv, {"start": c["start"], "end": c["end"], "inserted_bases": replacement}
    )
    target = jcv.annotations[c["feature_index"]]
    if any(a is target for a in dropped):
        return

    j = origin.index(c["feature_index"])
    exported = reread(write_genbank([new]))
    feat = exported.features[j]
    declared = feat.qualifiers["translation"][0]
    actual = str(Seq(str(feat.extract(exported.seq))).translate(table=1, to_stop=True))
    assert declared == actual, (
        f"{feat.location} still declares the pre-edit protein after "
        f"[{c['start']}, {c['end']}) was replaced inside it"
    )
