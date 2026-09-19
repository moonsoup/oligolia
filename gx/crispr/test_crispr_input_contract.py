"""CR-13, CR-14 — what the two endpoints refuse, and what a circle is.

CR-13 is the #87 family: one capability, two endpoints, one verdict on an input
that is not a sequence. CR-14's oracle is the same PAM-grammar scan run on the
circle instead of on the cut-open string.
"""

from __future__ import annotations

import pytest
from fastapi import HTTPException

from crispr_oracle import requirement, scan, scan_circular

CIRCULAR_INPUTS = ["puc19_L09137.2", "NC_001699.1", "AF071878.1"]


# ── CR-13: the two endpoints must refuse the same input ──────────────────────

# A space or a newline is layout, not a character: /crispr/design strips both
# before validating, and so should anything else. These four are characters that
# are not nucleotides, which is a different thing.
NOT_SEQUENCES = [
    pytest.param("ACGTACGTACGTACGTACG!", id="punctuation"),
    pytest.param("ACGTACGTACGTACGTACGU", id="rna-u"),
    pytest.param("ACGTACGTACGT12ACGTAC", id="digits"),
    pytest.param("ACGTACGTACGTACGT\tACG", id="tab"),
]


@pytest.mark.parametrize("guide", NOT_SEQUENCES)
def test_cr_13_design_refuses_a_non_dna_target(puc19, design, guide):
    """The reference behaviour: /crispr/design already refuses these."""
    target = puc19[:1000] + guide
    with pytest.raises(HTTPException) as excinfo:
        design(target, "SpCas9")
    assert 400 <= excinfo.value.status_code < 500


@pytest.mark.parametrize("guide", NOT_SEQUENCES)
def test_cr_13_score_guide_refuses_the_same_alphabet(guide):
    from backend.routers.crispr import score_guide

    params = requirement("CR-13")["parameters"]
    assert "17..24" in params["score_guide_length_guard"]
    assert 17 <= len(guide) <= 24, "length guard must not be what rejects it"

    try:
        got = score_guide(guide)
    except HTTPException as exc:
        assert 400 <= exc.status_code < 500
        return
    pytest.fail(
        "/crispr/score_guide accepted a guide that /crispr/design refuses and "
        f"answered with a verdict: {got['recommendation']!r}, on_target_score "
        f"{got['on_target_score']}, gc_content {got['gc_content']}"
    )


# ── CR-14: origin-spanning sites on a circular molecule ──────────────────────

@pytest.mark.parametrize("input_name", CIRCULAR_INPUTS)
@pytest.mark.parametrize("cas", ["SpCas9", "AsCas12a"])
def test_cr_14_origin_spanning_sites_are_found(corpus, design, input_name, cas):
    expected = requirement("CR-14")["parameters"]["expected"][input_name]
    target = corpus[input_name]

    linear = len(scan(cas, target))
    circular = scan_circular(cas, target)
    assert linear == expected["linear"][cas], "oracle disagrees with the register"
    assert circular == expected["circular"][cas], "oracle disagrees with the register"

    got = design(target, cas, max_guides=1)
    missed = circular - got.total_candidates
    assert got.total_candidates == circular, (
        f"{input_name} is a circular molecule; {cas} finds {circular} sites on the "
        f"circle but the app reports {got.total_candidates}, missing the {missed} "
        "site(s) that run through the origin"
    )
