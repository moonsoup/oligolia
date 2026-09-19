"""CO-2, CO-3, CO-4, CO-5, CO-6, CO-7 — what the app writes back out.

Export is the other half of the conversion: a 0-based half-open interval (s, e)
has to come back as ``(s+1)..e``, and every operator, partial boundary and
remote reference in the original location has to survive. The oracle is
``Bio.SeqIO`` re-reading the app's own output and comparing against
``Bio.SeqIO`` reading the pinned original — never the app's output against
itself.

Subject: ``backend/formats/genbank.py:write_genbank``, reached from the GUI's
"Save all sequences as GenBank…" (``gui/main_window.py:210``) and from
``POST /files/download/genbank``.
"""

from __future__ import annotations

import pytest
from Bio.SeqFeature import AfterPosition, BeforePosition

from conftest import (
    JCV,
    bio_extract,
    bio_parts,
    bio_profile,
    bio_record,
    bio_strand,
    feature_lines,
    genbank_with,
    requirement,
    reread,
)

CO_2 = requirement("CO-2")
CO_3 = requirement("CO-3")
CO_4 = requirement("CO-4")
CO_5 = requirement("CO-5")
CO_6 = requirement("CO-6")
CO_7 = requirement("CO-7")

INPUTS = [p.rsplit("/", 1)[-1] for p in CO_2["parameters"]["inputs"]]


def _round_trip(name):
    """(original record, record re-read from the subject's export)."""
    from backend.formats.genbank import read_genbank, write_genbank

    path = JCV.parent / name
    text = path.read_text()
    exported = write_genbank(read_genbank(text))
    return bio_record(path), reread(exported), exported


@pytest.mark.parametrize("name", INPUTS)
def test_co_2_round_trip_preserves_every_interval_and_strand(name):
    before, after, _ = _round_trip(name)
    assert [(f.type, bio_parts(f), bio_strand(f)) for f in after.features] == [
        (f.type, bio_parts(f), bio_strand(f)) for f in before.features
    ]


@pytest.mark.parametrize("name", INPUTS)
def test_co_2_round_trip_preserves_every_feature_s_bases(name):
    before, after, _ = _round_trip(name)
    assert str(after.seq) == str(before.seq)
    mismatched = [
        (f.type, str(f.location))
        for f, g in zip(before.features, after.features)
        if bio_extract(f, before) != bio_extract(g, after)
    ]
    assert mismatched == [], f"{len(mismatched)} features changed bases"


@pytest.mark.parametrize("name", INPUTS)
def test_co_2_exported_locations_are_one_based_inclusive(name):
    """(s, e) in memory must be written as (s+1)..e, which is what re-reading
    the export and subtracting one again has to give back."""
    before, after, _ = _round_trip(name)
    for f, g in zip(before.features, after.features):
        for (bs, be), (as_, ae) in zip(bio_parts(f), bio_parts(g)):
            assert (as_ + 1, ae) == (bs + 1, be)


@pytest.mark.parametrize("case", CO_3["parameters"]["cases"],
                         ids=lambda c: c["input"].rsplit("/", 1)[-1])
def test_co_3_origin_spanning_join_survives_export(case):
    _, after, exported = _round_trip(case["input"].rsplit("/", 1)[-1])
    feat = after.features[case["feature_index"]]
    assert bio_parts(feat) == [tuple(p) for p in case["parts"]]
    key, loc = feature_lines(exported)[case["feature_index"]]
    assert (key, loc) == (case["type"], case["file_location"])


@pytest.mark.parametrize("case", CO_3["parameters"]["cases"],
                         ids=lambda c: c["input"].rsplit("/", 1)[-1])
def test_co_3_exported_record_is_still_circular(case):
    _, after, _ = _round_trip(case["input"].rsplit("/", 1)[-1])
    assert after.annotations.get("topology") == CO_3["parameters"]["topology_after_export"]


def test_co_4_complement_join_is_written_back_as_complement_join():
    p = CO_4["parameters"]
    before, after, exported = _round_trip("NC_001699.1.gb")
    key, loc = feature_lines(exported)[p["feature_index"]]
    assert loc == p["exported_location_expected"]
    feat = after.features[p["feature_index"]]
    assert bio_strand(feat) == p["strand"]
    assert bio_parts(feat) == [tuple(q) for q in p["in_memory_parts"]]
    assert bio_extract(feat, after) == bio_extract(before.features[p["feature_index"]], before)


@pytest.mark.parametrize("case", CO_5["parameters"]["cases"],
                         ids=lambda c: f"{c['input'].rsplit('/', 1)[-1]}-{c['feature_index']}-{c['type']}")
def test_co_5_partial_boundaries_survive_export(case):
    """`<108..1007` must not come back as `108..1007` — INSDC gives `<` its own
    meaning (the feature starts before the first sequenced base)."""
    before, after, exported = _round_trip(case["input"].rsplit("/", 1)[-1])
    original = before.features[case["feature_index"]].location
    exported_loc = after.features[case["feature_index"]].location
    was_fuzzy = (isinstance(original.parts[0].start, BeforePosition)
                 or isinstance(original.parts[-1].end, AfterPosition))
    assert was_fuzzy, "fixture drift: this feature is not partial in the pinned file"
    still_fuzzy = (isinstance(exported_loc.parts[0].start, BeforePosition)
                   or isinstance(exported_loc.parts[-1].end, AfterPosition))
    assert still_fuzzy, (
        f"{case['file_location']} exported as "
        f"{feature_lines(exported)[case['feature_index']][1]!r}"
    )


@pytest.mark.parametrize("case", CO_5["parameters"]["cases"],
                         ids=lambda c: f"{c['input'].rsplit('/', 1)[-1]}-{c['feature_index']}-{c['type']}")
def test_co_5_exported_location_string_still_carries_the_angle_bracket(case):
    _, _, exported = _round_trip(case["input"].rsplit("/", 1)[-1])
    _, loc = feature_lines(exported)[case["feature_index"]]
    assert "<" in loc or ">" in loc, f"{case['file_location']} -> {loc!r}"


def test_co_6_order_is_not_rewritten_as_join():
    """INSDC: join asserts the parts form one contiguous sequence; order does
    not. Rewriting one as the other changes the record's claim."""
    from backend.formats.genbank import read_genbank, write_genbank

    p = CO_6["parameters"]
    text = genbank_with([
        ("misc_feature", p["probe_location"], "ordered"),
        ("misc_feature", p["companion_control"], "joined"),
    ])
    assert reread(text).features[1].location.operator == "order"  # fixture sanity
    exported = write_genbank(read_genbank(text))
    assert reread(exported).features[1].location.operator == p["expected_exported_operator"], (
        f"{p['probe_location']} exported as {feature_lines(exported)[1][1]!r}"
    )


def test_co_6_a_real_join_still_exports_as_join():
    """The control: fixing order() must not be done by breaking join()."""
    from backend.formats.genbank import read_genbank, write_genbank

    p = CO_6["parameters"]
    text = genbank_with([
        ("misc_feature", p["probe_location"], "ordered"),
        ("misc_feature", p["companion_control"], "joined"),
    ])
    exported = write_genbank(read_genbank(text))
    assert reread(exported).features[2].location.operator == "join"


def test_co_7_remote_reference_is_not_silently_localised():
    """`J00194.1:100..202` names bases in another entry. Exporting it as a bare
    `100..202` hands this record's own bases to the feature instead."""
    from backend.formats.genbank import read_genbank, write_genbank

    p = CO_7["parameters"]
    text = genbank_with([("misc_feature", p["probe_location"], "remote")])
    assert reread(text).features[1].location.parts[1].ref == p["remote_accession"]

    try:
        exported = write_genbank(read_genbank(text))
    except Exception:
        return  # rejecting the feature is an acceptable outcome
    _, loc = feature_lines(exported)[1]
    assert loc != p["forbidden_export"], (
        f"{p['probe_location']} exported as {loc!r}: bases 100..202 of this "
        "record now stand in for bases of J00194.1"
    )
    assert p["remote_accession"] in loc


def test_co_7_the_local_part_of_the_remote_join_is_unharmed():
    """Whatever happens to the remote part, 201..210 is this record's own."""
    from backend.formats.genbank import read_genbank

    p = CO_7["parameters"]
    text = genbank_with([("misc_feature", p["probe_location"], "remote")])
    ann = read_genbank(text)[0].annotations[1]
    assert tuple(ann.parts[0]) == (200, 210)


def test_co_2_no_feature_is_added_or_dropped():
    """A round trip that quietly loses a feature would pass every zip() above."""
    for name in INPUTS:
        before, after, _ = _round_trip(name)
        assert len(after.features) == len(before.features), name
        assert len(bio_profile(after)) == len(bio_profile(before))
