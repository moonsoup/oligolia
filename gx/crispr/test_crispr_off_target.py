"""CR-11, CR-12 — the specificity number, and whether SpCas9-HF1 is a nuclease.

CR-11's oracle is the Hsu et al. 2013 / MIT definition re-typed from
gx/corpus/crispr/crispor-v3.1.py (`hitScoreM`, `calcHitScore`,
`calcMitGuideScore`) and checked against those bytes at import time. CR-12 needs
no external authority — the CRISPR panel's own Cas-type description is the claim
under test.
"""

from __future__ import annotations

import pytest

from crispr_oracle import HIT_SCORE_M, mit_guide_score, mit_hit_score, requirement


# ── CR-11: an MIT score that ignores where the mismatch is ───────────────────

def _carrier(protospacer: str) -> str:
    """One NGG site, in flanks that form no second site (see CR-11 parameters)."""
    return "TTTTTT" + protospacer + "AGG" + "TTTTTT"


def _mutate(guide: str, position_1based: int) -> str:
    bases = list(guide)
    i = position_1based - 1
    bases[i] = "A" if bases[i] != "A" else "C"
    return "".join(bases)


@pytest.mark.parametrize("position", [20, 1])
def test_cr_11_specificity_depends_on_where_the_mismatch_sits(position):
    from backend.crispr_offtarget import scan_off_targets

    params = requirement("CR-11")["parameters"]
    guide = params["guide"]
    key = f"mismatch_at_position_{position}"
    expected = params["expected"][key]

    off_target = _mutate(guide, position)
    assert sum(a != b for a, b in zip(guide, off_target)) == 1

    # the oracle, from the pinned definition
    hit = mit_hit_score(guide, off_target)
    assert hit == pytest.approx(expected["hit_score"], abs=0.1)
    assert mit_guide_score(hit) == pytest.approx(expected["specificity_score"], abs=0.1)
    assert HIT_SCORE_M[19] == 0.583 and HIT_SCORE_M[0] == 0.0

    result = scan_off_targets(
        guide,
        [_carrier(guide), _carrier(off_target)],
        cas_family=params["cas_family"],
        max_mismatches=params["max_mismatches"],
    )
    assert result.summary["1"] == 1, f"expected exactly one 1-mismatch hit, got {result.summary}"
    assert result.specificity_score == pytest.approx(expected["specificity_score"], abs=0.1), (
        f"a single mismatch at protospacer position {position} (Hsu weight "
        f"{HIT_SCORE_M[position - 1]}) should give Sguide "
        f"{expected['specificity_score']}, the app reports {result.specificity_score}"
    )


def test_cr_11_two_mismatch_positions_do_not_give_the_same_specificity():
    """The single sharpest statement of CR-11, in one assertion."""
    from backend.crispr_offtarget import scan_off_targets

    params = requirement("CR-11")["parameters"]
    guide = params["guide"]
    seed, distal = _mutate(guide, 20), _mutate(guide, 1)

    assert mit_guide_score(mit_hit_score(guide, seed)) != pytest.approx(
        mit_guide_score(mit_hit_score(guide, distal)), abs=0.1
    ), "oracle: the pinned definition must separate these two"

    a = scan_off_targets(guide, [_carrier(guide), _carrier(seed)], max_mismatches=3)
    b = scan_off_targets(guide, [_carrier(guide), _carrier(distal)], max_mismatches=3)
    assert a.specificity_score != b.specificity_score, (
        "a PAM-proximal seed mismatch and a PAM-distal one both score "
        f"{a.specificity_score}; the app's specificity is a function of the "
        "mismatch-bucket counts alone and never reads the positions"
    )


# ── CR-12: is SpCas9-HF1 anything at all? ────────────────────────────────────

@pytest.mark.parametrize("input_name", ["puc19_L09137.2", "NG_007400.1"])
def test_cr_12_sp_cas9_hf1_is_not_an_alias_of_sp_cas9(corpus, design, input_name):
    from gui.panels.crispr_panel import CAS_INFO

    assert "fewer off-targets" in CAS_INFO["SpCas9-HF1"]

    target = corpus[input_name]
    kw = {"max_guides": 10, "check_off_targets": True}
    plain = design(target, "SpCas9", **kw)
    hifi = design(target, "SpCas9-HF1", **kw)

    same = (
        plain.total_candidates == hifi.total_candidates
        and [g.model_dump() for g in plain.guides] == [g.model_dump() for g in hifi.guides]
    )
    assert not same, (
        "SpCas9-HF1 returns a field-for-field identical response to SpCas9 "
        f"({plain.total_candidates} candidates, {len(plain.guides)} guides, "
        f"identical specificity scores {sorted({g.specificity_score for g in plain.guides})}) "
        "while the panel describes it as 'high fidelity · fewer off-targets'"
    )
