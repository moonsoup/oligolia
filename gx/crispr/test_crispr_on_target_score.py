"""CR-8, CR-9, CR-10 — the number in the "On-target" column, and what it picks.

The oracle for CR-8 is the pinned Doench et al. 2014 Rule Set 1 definition,
re-typed from gx/corpus/crispr/crisporEffScores-v3.1.py and checked against
those bytes at import time. CR-9 and CR-10 need no external authority: they
hold the app to its own single score field and to its own candidate count.
"""

from __future__ import annotations

import pytest

from crispr_oracle import (
    CANONICAL_LENGTH,
    doench_rule_set1,
    mer30_for,
    overlaps,
    requirement,
    scan,
    spearman,
)


# ── CR-8: the score named after Rule Set 1 ───────────────────────────────────

@pytest.fixture(scope="module")
def rs1_vs_app(puc19):
    """(protospacer, Rule Set 1, app score) for every scorable pUC19 SpCas9 site.

    Rule Set 1 is evaluated on each site's OWN real 30mer, cut out of pUC19 by
    index; sites too close to an end for the 4 nt of 5' context or 3 nt of 3'
    context are skipped rather than padded, because the model is undefined there.
    """
    from backend.routers.crispr import score_guide

    rows = []
    for strand, pos, proto, _pam in scan("SpCas9", puc19):
        mer30 = mer30_for(puc19, pos, strand)
        if mer30 is None:
            continue
        rows.append((proto, doench_rule_set1(mer30), score_guide(proto)["on_target_score"]))
    return rows


def test_cr_8_reported_score_reproduces_rule_set_1_on_the_real_30mer(rs1_vs_app):
    """Every pUC19 SpCas9 site, scored by the app and by the pinned definition."""
    params = requirement("CR-8")["parameters"]
    assert len(rs1_vs_app) == params["scorable_sites"]
    assert params["rs1_30mer"].startswith("4 nt 5' context")

    off = [(proto, rs1, app) for proto, rs1, app in rs1_vs_app if abs(rs1 - app) > 0.05]
    worst = max(rs1_vs_app, key=lambda r: abs(r[1] - r[2]))
    assert not off, (
        f"{len(off)} of {len(rs1_vs_app)} sites disagree with pinned Rule Set 1 by "
        f"more than 0.05; worst is {worst[0]} — Rule Set 1 {worst[1]:.3f}, app "
        f"{worst[2]:.3f}. Rule Set 1 is defined on the 30mer (PAM and both flanks "
        "carry weights); the app scores the 20 nt protospacer alone."
    )


def test_cr_8_app_ordering_agrees_with_rule_set_1(rs1_vs_app):
    params = requirement("CR-8")["parameters"]
    assert len(rs1_vs_app) == params["scorable_sites"]

    rho = spearman([r[1] for r in rs1_vs_app], [r[2] for r in rs1_vs_app])
    assert rho >= params["min_spearman_rho"], (
        f"Spearman rho between the app's on-target score and pinned Rule Set 1 "
        f"over {len(rs1_vs_app)} real pUC19 sites is {rho:.4f}"
    )


def test_cr_8_max_guides_would_select_the_same_guides(rs1_vs_app):
    params = requirement("CR-8")["parameters"]
    top_rs1 = {r[0] for r in sorted(rs1_vs_app, key=lambda r: -r[1])[:10]}
    top_app = {r[0] for r in sorted(rs1_vs_app, key=lambda r: -r[2])[:10]}
    shared = len(top_rs1 & top_app)

    best_rs1 = max(rs1_vs_app, key=lambda r: r[1])
    assert best_rs1[0] == params["rs1_best_guide"]
    assert best_rs1[1] == pytest.approx(params["rs1_best_score"], abs=1e-4)

    assert shared >= params["min_top10_overlap"], (
        f"only {shared} of the app's top 10 guides are in Rule Set 1's top 10; "
        f"Rule Set 1's best guide is {best_rs1[0]} ({best_rs1[1]:.3f}), the app "
        f"ranks it {sorted(rs1_vs_app, key=lambda r: -r[2]).index(best_rs1) + 1}"
    )


# ── CR-9: one field, one scale ───────────────────────────────────────────────

@pytest.mark.parametrize("cas", ["SpCas9", "SpCas9-HF1", "AsCas12a", "LwaCas13a"])
def test_cr_9_design_and_score_guide_agree_on_the_same_guide(puc19, design, cas):
    from backend.routers.crispr import score_guide

    resp = design(puc19, cas, max_guides=10)
    mismatches = []
    for g in resp.guides:
        theirs = score_guide(g.sequence)["on_target_score"]
        if g.on_target_score != theirs:
            mismatches.append((g.sequence, g.on_target_score, theirs))
    assert not mismatches, (
        f"{cas}: /crispr/design and /crispr/score_guide disagree on "
        f"{len(mismatches)} of {len(resp.guides)} guides, e.g. {mismatches[0]}"
    )


# ── CR-10: does max_guides select, or collapse? ──────────────────────────────

@pytest.mark.parametrize(
    ("input_name", "cas"),
    [("NG_007400.1", "LwaCas13a"), ("puc19_L09137.2", "LwaCas13a"),
     ("NG_007400.1", "AsCas12a"), ("NG_007400.1", "SpCas9")],
)
def test_cr_10_shown_guides_are_not_one_locus_in_n_shifted_copies(
    corpus, design, input_name, cas
):
    params = requirement("CR-10")["parameters"]
    target = corpus[input_name]
    resp = design(target, cas, max_guides=params["max_guides"])
    positions = [g.position for g in resp.guides]
    L = CANONICAL_LENGTH[cas]

    non_overlapping = sum(
        1
        for i, a in enumerate(positions)
        for b in positions[i + 1:]
        if not overlaps(a, b, L)
    )
    assert non_overlapping >= params["min_non_overlapping_pairs"], (
        f"{cas} on {input_name}: {resp.total_candidates} candidates searched, "
        f"{len(positions)} shown, and every pair overlaps — positions {positions} "
        f"span {max(positions) - min(positions)} nt of a {len(target)} nt molecule"
    )
