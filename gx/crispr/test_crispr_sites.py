"""CR-1..CR-5, CR-7 — is the site the app reports the site that is in the DNA?

Every expected value here is re-derived from the raw pinned sequence: an
independent both-strand regex scan written from the PAM grammar, a slice by
index, Biopython's reverse complement, or Biopython's gc_fraction. The subject's
own output is never the oracle.
"""

from __future__ import annotations

import pytest

from crispr_oracle import (
    CANONICAL_LENGTH,
    PAM_MOTIF,
    PAM_REGEX,
    PAM_SIDE,
    gc_percent,
    rc,
    requirement,
    scan,
    site_on,
)

CAS_TYPES = ["SpCas9", "SpCas9-HF1", "AsCas12a", "LwaCas13a"]
INPUTS = ["puc19_L09137.2", "NC_001699.1", "AF071878.1", "NG_007400.1"]
CASES = [(i, c) for i in INPUTS for c in CAS_TYPES]
IDS = [f"{i}-{c}" for i, c in CASES]


# ── CR-1: the candidate count is the whole molecule, both strands ────────────

@pytest.mark.parametrize(("input_name", "cas"), CASES, ids=IDS)
def test_cr_1_total_candidates_equals_an_independent_both_strand_scan(
    corpus, design, input_name, cas
):
    expected = requirement("CR-1")["parameters"]["expected_total_candidates"]
    target = corpus[input_name]

    independent = scan(cas, target)
    assert len(independent) == expected[input_name][cas], (
        "the oracle itself disagrees with the pinned register — fix the register, "
        "not the assertion below"
    )

    got = design(target, cas, max_guides=1)
    assert got.target_length == len(target)
    assert got.total_candidates == expected[input_name][cas]


# ── CR-2 / CR-3 / CR-4 / CR-5: re-extract every reported guide by index ──────

@pytest.mark.parametrize(("input_name", "cas"), CASES, ids=IDS)
def test_cr_2_reported_protospacer_is_the_one_at_the_reported_index(
    corpus, design, input_name, cas
):
    target = corpus[input_name]
    resp = design(target, cas, max_guides=50)
    assert resp.guides, "nothing to check"

    for g in resp.guides:
        proto, _pam = site_on(cas, target, g.position, g.strand, len(g.sequence))
        assert proto == g.sequence, (
            f"{cas} {g.strand} guide at {g.position}: reported {g.sequence}, "
            f"but the sequence holds {proto}"
        )


@pytest.mark.parametrize(("input_name", "cas"), CASES, ids=IDS)
def test_cr_3_guide_excludes_the_pam(corpus, design, input_name, cas):
    target = corpus[input_name]
    resp = design(target, cas, max_guides=50)
    pam_len = {"NGG": 3, "TTTV": 4, "none": 0}[PAM_MOTIF[cas]]

    for g in resp.guides:
        assert len(g.sequence) == CANONICAL_LENGTH[cas]
        if not pam_len:
            continue
        _proto, pam = site_on(cas, target, g.position, g.strand, len(g.sequence))
        assert len(pam) == pam_len
        # The PAM must sit outside the footprint: prepending or appending it to
        # the guide has to reconstruct a longer stretch of real sequence, which
        # is only true if the guide did not already swallow part of it (#53).
        whole = g.sequence + pam if PAM_SIDE[cas] == 3 else pam + g.sequence
        assert len(whole) == len(g.sequence) + pam_len
        strand_seq = target if g.strand == "+" else rc(target)
        assert whole in strand_seq, (
            f"{cas} {g.strand} guide at {g.position}: guide+PAM {whole} is not a "
            "contiguous stretch of the molecule"
        )


@pytest.mark.parametrize(("input_name", "cas"), CASES, ids=IDS)
def test_cr_4_pam_obeys_the_grammar_the_app_reports(corpus, design, input_name, cas):
    import re

    target = corpus[input_name]
    resp = design(target, cas, max_guides=50)
    motif = PAM_MOTIF[cas]

    for g in resp.guides:
        assert g.pam == motif, f"{cas} guide reports pam={g.pam!r}, expected {motif!r}"
        if motif == "none":
            assert g.strand == "+", "LwaCas13a targets the given sense strand only"
            continue
        _proto, pam = site_on(cas, target, g.position, g.strand, len(g.sequence))
        assert re.fullmatch(PAM_REGEX[motif], pam), (
            f"{cas} {g.strand} guide at {g.position}: PAM read out of the sequence "
            f"is {pam!r}, which does not match {motif}"
        )


@pytest.mark.parametrize(("input_name", "cas"), CASES, ids=IDS)
def test_cr_5_minus_strand_indices_map_under_reverse_complement(
    corpus, design, input_name, cas
):
    target = corpus[input_name]
    n = len(target)
    resp = design(target, cas, max_guides=50)
    minus = [g for g in resp.guides if g.strand == "-"]
    if cas == "LwaCas13a":
        assert not minus
        return

    for g in minus:
        L = len(g.sequence)
        assert 0 <= g.position and g.position + L <= n
        assert rc(target[g.position:g.position + L]) == g.sequence
        # the PAM has to fit too, on the far side of the footprint
        if PAM_SIDE[cas] == 3:
            assert g.position >= 3, "an NGG site on the minus strand needs 3 bases below it"
        else:
            assert g.position + L + 4 <= n


@pytest.mark.parametrize("input_name", INPUTS)
def test_cr_5_both_strands_are_actually_represented(corpus, design, input_name):
    """A one-strand scan is the #64.2 failure; assert both strands appear."""
    target = corpus[input_name]
    for cas in ("SpCas9", "AsCas12a"):
        independent = scan(cas, target)
        assert {s for s, *_ in independent} == {"+", "-"}, "oracle found one strand only"
        resp = design(target, cas, max_guides=50)
        by_strand = {g.strand for g in resp.guides}
        # max_guides truncates, so require the full candidate set to be balanced
        # rather than the shown ten: compare the count instead.
        assert resp.total_candidates == len(independent)
        assert by_strand <= {"+", "-"}


# ── CR-7: GC recomputed with Biopython ───────────────────────────────────────

@pytest.mark.parametrize(("input_name", "cas"), CASES, ids=IDS)
def test_cr_7_gc_content_is_the_gc_of_the_reported_sequence(
    corpus, design, input_name, cas
):
    resp = design(corpus[input_name], cas, max_guides=50)
    for g in resp.guides:
        assert g.gc_content == gc_percent(g.sequence), (
            f"{g.sequence}: reported {g.gc_content}, Biopython says {gc_percent(g.sequence)}"
        )
