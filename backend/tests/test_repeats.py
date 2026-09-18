"""#65: find_repeats reported fragments instead of repeats, and dropped real ones.

Three defects:
  1. a single (CA)40 microsatellite came back as ~50 overlapping shifted fragments
     of itself, because the dedup key was (start, unit) and every offset into the
     run is a different start;
  2. a genuine (GGT)6 elsewhere in the same sequence was not reported at all,
     because `unique[:50]` truncated the list BEFORE it was sorted by position;
  3. the docstring and the `repeat_type` field advertised inverted repeats, which
     were never implemented.
"""

from __future__ import annotations

from fastapi.testclient import TestClient


def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _repeats(client: TestClient, seq: str, **params) -> list[dict]:
    r = client.post("/analysis/find_repeats", params={"sequence": seq, **params})
    assert r.status_code == 200, r.text
    return r.json()


def test_one_microsatellite_is_one_repeat(client: TestClient) -> None:
    """#65.1: (CA)40 was returned as ~50 shifted fragments of itself."""
    seq = "GGGATC" + "CA" * 40 + "GGGATC"
    tandem = [x for x in _repeats(client, seq) if x["repeat_type"] == "tandem"]

    ca = [x for x in tandem if set(x["unit"]) <= {"C", "A"}]
    assert len(ca) == 1, f"expected one CA repeat, got {len(ca)}: {ca[:5]}"
    hit = ca[0]
    assert hit["unit"] == "CA", hit
    assert hit["start"] == 6, hit
    assert hit["end"] == 6 + 80, hit
    assert hit["copies"] == 40, hit


def test_a_second_real_repeat_is_not_dropped(client: TestClient) -> None:
    """#65.2: the 50-item cap was applied before sorting, losing later repeats."""
    seq = "GGGATC" + "CA" * 40 + "TTTTAAAACCCC" + "GGT" * 6 + "GGGATC"
    tandem = [x for x in _repeats(client, seq) if x["repeat_type"] == "tandem"]

    units = {x["unit"] for x in tandem}
    assert "CA" in units, units
    assert "GGT" in units, f"the (GGT)6 repeat was dropped: {units}"

    ggt = [x for x in tandem if x["unit"] == "GGT"][0]
    assert ggt["copies"] == 6, ggt


def test_the_unit_reported_is_the_primitive_one(client: TestClient) -> None:
    """(CACA)20 and (CA)40 are the same repeat; report the primitive unit."""
    seq = "GGGATC" + "CA" * 40 + "GGGATC"
    tandem = [x for x in _repeats(client, seq) if x["repeat_type"] == "tandem"]
    ca = [x for x in tandem if set(x["unit"]) <= {"C", "A"}]
    assert [x["unit"] for x in ca] == ["CA"], [x["unit"] for x in ca]


def test_no_repeat_is_contained_in_another_and_none_share_a_unit(client: TestClient) -> None:
    """The property that makes 'fragments' impossible, stated correctly.

    NOT "no two repeats overlap": adjacent runs can legitimately share a base.
    In "…CACA" + "ATAT…"*8 + "GGT"*8 the AT run ends at 82 and a real TGG run
    starts at 81, because the AT run's final T begins the TGG periodicity. That
    is a true repeat, not a fragment.

    What must never happen is the same locus reported many times: one repeat
    inside another, or two runs of the same unit overlapping.
    """
    seq = "GGGATC" + "CA" * 30 + "ATATATATATATATAT" + "GGT" * 8 + "TTTTT"
    tandem = [x for x in _repeats(client, seq) if x["repeat_type"] == "tandem"]

    for a in tandem:
        for b in tandem:
            if a is b:
                continue
            contained = b["start"] <= a["start"] and a["end"] <= b["end"]
            assert not contained, f"{a} is contained in {b}"
            if a["unit"] == b["unit"]:
                assert a["end"] <= b["start"] or b["end"] <= a["start"], (a, b)


def test_results_are_sorted_by_position(client: TestClient) -> None:
    seq = "GGGATC" + "CA" * 20 + "TTTT" + "GGT" * 6 + "AAAA" + "AT" * 15
    out = _repeats(client, seq)
    assert [x["start"] for x in out] == sorted(x["start"] for x in out)


def test_an_inverted_repeat_is_found(client: TestClient) -> None:
    """#65.3: advertised by the docstring and the repeat_type field, never implemented.

    The spacer is deliberately A-only. A GC spacer is self-complementary, so the
    arms genuinely extend into it and the implementation correctly reports a
    LONGER inverted repeat than the one planted — which is right, and made an
    exact-unit assertion the wrong thing to write.
    """
    arm = "AAGGCCTTAG"
    spacer = "AAAAAAAA"
    seq = "GGGG" + arm + spacer + _rc(arm) + "GGGG"
    inverted = [x for x in _repeats(client, seq) if x["repeat_type"] == "inverted"]

    assert inverted, "no inverted repeat found in a sequence built to contain one"
    hit = max(inverted, key=lambda x: x["length"])
    assert hit["unit"] == arm, hit
    assert hit["start"] == 4, hit
    assert hit["end"] == 4 + len(arm) + len(spacer) + len(arm), hit


def test_a_reported_inverted_repeat_is_always_really_one(client: TestClient) -> None:
    """Whatever arm it settles on, the arms must be reverse complements.

    This is the assertion that holds even when the arms extend into a
    self-complementary spacer, which is the case that broke a naive test.
    """
    arm = "AAGGCCTTAG"
    for spacer in ("AAAAAAAA", "GCGCGCGC", "ACGTACGT", ""):
        seq = "GGGG" + arm + spacer + _rc(arm) + "GGGG"
        for hit in _repeats(client, seq):
            if hit["repeat_type"] != "inverted":
                continue
            unit = hit["unit"]
            left = seq[hit["start"]:hit["start"] + len(unit)]
            right = seq[hit["end"] - len(unit):hit["end"]]
            assert left == unit, (spacer, hit)
            assert left == _rc(right), (spacer, hit, left, right)


def test_a_sequence_with_no_inverted_repeat_reports_none(client: TestClient) -> None:
    seq = "AAAAAAAAAACCCCCCCCCCAAAAAAAAAACCCCCCCCCC"
    inverted = [x for x in _repeats(client, seq) if x["repeat_type"] == "inverted"]
    for hit in inverted:
        # Whatever is reported must actually BE an inverted repeat.
        left = seq[hit["start"]:hit["start"] + len(hit["unit"])]
        right = seq[hit["end"] - len(hit["unit"]):hit["end"]]
        assert left == _rc(right), hit


def test_every_reported_repeat_re_derives_from_the_sequence(client: TestClient) -> None:
    """The re-derivation invariant, as for CRISPR guides."""
    seq = "GGGATC" + "CA" * 25 + "TTTT" + "GGT" * 7 + "AAGGCCTTAG" + "GC" * 4 + "CTAAGGCCTT"
    for hit in _repeats(client, seq):
        span = seq[hit["start"]:hit["end"]]
        assert len(span) == hit["length"], hit
        if hit["repeat_type"] == "tandem":
            unit = hit["unit"]
            assert span.startswith(unit), hit
            assert span == unit * int(hit["copies"]), hit
        else:
            unit = hit["unit"]
            assert span.startswith(unit), hit
            assert span.endswith(_rc(unit)), hit


def test_short_sequences_do_not_crash(client: TestClient) -> None:
    for seq in ("", "A", "AC", "ACGT"):
        r = client.post("/analysis/find_repeats", params={"sequence": seq})
        assert r.status_code in (200, 400), (seq, r.text)
