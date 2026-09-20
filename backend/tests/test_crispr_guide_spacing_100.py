"""#100: max_guides must select among the candidates, not slice one locus.

The Cas13 on-target score takes two values and the GC tie-break is exactly 0.0
for every 50%-GC window, so on a long molecule the whole candidate set ties and
a stable sort used to return the first `max_guides` windows of the target —
ten 22-mers 14 nt apart, while the panel reported "10 guides shown (of 24,523
candidates)".
"""

from fastapi.testclient import TestClient

from ._puc19 import PUC19

# Every 22-mer of an "AG" repeat is 11 A + 11 G, i.e. exactly 50% GC, so both
# the score and the GC tie-break are constant over the entire candidate set.
# This is the synthetic form of the COL1A1 case in #100.
TIED_TARGET = "AG" * 1000


def _positions(data: dict) -> list[int]:
    return [g["position"] for g in data["guides"]]


def _non_overlapping_pairs(positions: list[int], length: int) -> int:
    return sum(
        1
        for i, a in enumerate(positions)
        for b in positions[i + 1:]
        if abs(a - b) >= length
    )


def test_cas13_shown_guides_are_not_one_locus_in_shifted_copies(
    client: TestClient,
) -> None:
    r = client.post("/crispr/design", json={
        "target_sequence": TIED_TARGET, "cas_type": "LwaCas13a", "max_guides": 10,
    })
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["total_candidates"] == len(TIED_TARGET) - 22 + 1
    positions = _positions(data)
    assert len(positions) == 10
    assert _non_overlapping_pairs(positions, 22) >= 1, (
        f"every shown pair overlaps: positions {positions} span "
        f"{max(positions) - min(positions)} nt of a {len(TIED_TARGET)} nt target"
    )
    # With 1,979 tied candidates there is no reason for ANY shown pair to overlap.
    assert _non_overlapping_pairs(positions, 22) == 45, positions


def test_cas9_shown_guides_do_not_overlap_on_a_real_plasmid(
    client: TestClient,
) -> None:
    r = client.post("/crispr/design", json={
        "target_sequence": PUC19, "cas_type": "SpCas9", "max_guides": 10,
    })
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["total_candidates"] > 100
    positions = _positions(data)
    assert len(positions) == 10
    assert _non_overlapping_pairs(positions, 20) >= 1, positions


def test_spacing_does_not_shrink_the_shown_count_on_a_short_target(
    client: TestClient,
) -> None:
    """A 30 nt target has 9 Cas13 windows, all mutually overlapping.

    The spacing rule must not turn "5 guides shown" into "1 guide shown": the
    remaining slots are filled from the rest of the ranking.
    """
    target = "AG" * 15
    r = client.post("/crispr/design", json={
        "target_sequence": target, "cas_type": "LwaCas13a", "max_guides": 5,
    })
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["total_candidates"] == 9
    assert len(data["guides"]) == 5
    assert len(set(_positions(data))) == 5, "the same window must not be shown twice"


def test_min_guide_spacing_zero_restores_the_plain_slice(
    client: TestClient,
) -> None:
    """The old behaviour stays reachable, explicitly, via the new parameter."""
    r = client.post("/crispr/design", json={
        "target_sequence": TIED_TARGET, "cas_type": "LwaCas13a",
        "max_guides": 10, "min_guide_spacing": 0,
    })
    assert r.status_code == 200, r.text
    assert _positions(r.json()) == list(range(10))


def test_shown_guides_stay_in_rank_order(client: TestClient) -> None:
    """Spreading the set out must not scramble the on-target ranking."""
    r = client.post("/crispr/design", json={
        "target_sequence": PUC19, "cas_type": "SpCas9", "max_guides": 10,
    })
    assert r.status_code == 200, r.text
    scores = [g["on_target_score"] for g in r.json()["guides"]]
    assert scores == sorted(scores, reverse=True), scores
