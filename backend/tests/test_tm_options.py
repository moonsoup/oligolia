"""Tm method and buffer are selectable, and every pair says which it used (#81).

Split out of #79: the *mechanical* half. Before this, `TM_CONDITIONS` was pinned
inside `backend/routers/primers.py` with no way for a caller to ask for anything
else, and `_tm_wallace` — described by #50 as "kept as an option" — had no request
parameter reaching it, so it was a retained function rather than an option.

What this file pins:

  * the DEFAULT is unchanged. `DEFAULT_PAIRS` is recorded from `main` at d1ebad1's
    successor (487da5d) BEFORE this change, so a request that sets none of the new
    fields has to reproduce it exactly. Choosing a different default buffer is #79's
    open scientific question and is explicitly out of scope here.
  * `tm_method="wallace"` really reaches `_tm_wallace`.
  * a PCR-like buffer really reaches `Tm_NN`, against values RECORDED from
    Biopython 1.85 on 2026-09-18 — not computed live in this test.
  * buffers `Tm_NN` cannot evaluate are refused with a 4xx rather than crashing
    with `ValueError: math domain error` (a 500).

The recorded numbers below are the whole point of #79's evidence: about 10 degC
separates the pinned bench buffer from a buffer someone actually runs.
"""

from __future__ import annotations

import pytest

from backend.routers.primers import TM_CONDITIONS, _tm_nearest_neighbor, _tm_wallace

#: The M13 reverse primer, and its Tm RECORDED from Biopython 1.85 (DNA_NN3,
#: saltcorr=5) on 2026-09-18 under two buffers. Nothing here calls Tm_NN.
M13_REVERSE = "AGCGGATAACAATTTCACACAGG"
M13_TM_DEFAULT_BUFFER = 53.7   # Na 50 mM, Mg 0, primer 25 nM, template 25 nM
M13_TM_PCR_BUFFER = 64.0       # Na 50 mM, Mg 1.5 mM, dNTPs 0.2 mM, primer 250 nM, template 0

#: A template known to yield pairs. Used for the default-unchanged pin.
TEMPLATE = (
    "ATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACGTGGCATCACGTTT"
    "CGGCCATGGACGTAGGCATCACGTGGCATCACGTGGCCGGCGGCGGCGCTGCTGCAAGCTT"
    "ACGTGGCATCACGATGGCCTGTGGGCATTTGGCCAATTTAGGCCATGGACGTGGCATCACG"
)

#: The same template with the M13 reverse primer spliced on the front, so the
#: endpoint can be asked for a primer whose Tm is independently recorded above.
M13_TEMPLATE = M13_REVERSE + TEMPLATE

DEFAULT_REQUEST = {
    "template": TEMPLATE,
    "tm_min": 45.0,
    "tm_max": 75.0,
    "product_min": 80,
    "product_max": 200,
    "max_pairs": 5,
}

#: Recorded from `main` BEFORE #81, by POSTing DEFAULT_REQUEST. This is the
#: "defaults do not move" pin: #81's out-of-scope line says the default buffer
#: stays exactly TM_CONDITIONS, and this is what that means in observable terms.
DEFAULT_PAIRS = [
    {"forward": "CCTGTGGGCATTTGGCCA", "forward_tm": 54.8,
     "reverse": "CCAAATGCCCACAGGCCA", "reverse_tm": 54.8,
     "product_size": 150, "penalty": 0.0},
    {"forward": "TGTGGGCATTTGGCCAAT", "forward_tm": 51.9,
     "reverse": "ATGCCACGTAAGCTTGCA", "reverse_tm": 51.9,
     "product_size": 125, "penalty": 0.0},
    {"forward": "AGGCCATGGACGTGGCAT", "forward_tm": 55.8,
     "reverse": "ATGCCCACAGGCCATCGT", "reverse_tm": 55.8,
     "product_size": 124, "penalty": 0.0},
    {"forward": "GCCATGGACGTGGCATCA", "forward_tm": 54.7,
     "reverse": "AATGCCCACAGGCCATCG", "reverse_tm": 54.7,
     "product_size": 123, "penalty": 0.0},
    {"forward": "GTGGCATCACGTGGCATC", "forward_tm": 53.6,
     "reverse": "CCCACAGGCCATCGTGAT", "reverse_tm": 53.6,
     "product_size": 110, "penalty": 0.0},
]

#: Narrow GC window that lets the M13 reverse primer (43.5% GC) through while
#: keeping the GC-rich rest of the template out, so it reaches the returned set.
M13_REQUEST = {
    "template": M13_TEMPLATE,
    "tm_min": 45.0,
    "tm_max": 75.0,
    "gc_min": 40.0,
    "gc_max": 48.0,
    "product_min": 80,
    "product_max": 200,
    "primer_len_min": 23,
    "primer_len_max": 23,
    "max_pairs": 20,
}

PCR_BUFFER = {"mg_mm": 1.5, "dntp_mm": 0.2, "primer_nm": 250, "template_nm": 0}


def _post(client, **overrides):
    body = dict(DEFAULT_REQUEST)
    body.update(overrides)
    return client.post("/primers/design", json=body)


def _summarise(pairs: list[dict]) -> list[dict]:
    return [
        {"forward": p["forward"]["sequence"], "forward_tm": p["forward"]["tm"],
         "reverse": p["reverse"]["sequence"], "reverse_tm": p["reverse"]["tm"],
         "product_size": p["product_size"], "penalty": p["penalty"]}
        for p in pairs
    ]


def _tm_by_sequence(pairs: list[dict]) -> dict[str, float]:
    return {p[role]["sequence"]: p[role]["tm"]
            for p in pairs for role in ("forward", "reverse")}


# --- acceptance 1: a default request returns the same pairs and Tms as main ----

def test_a_default_request_is_byte_for_byte_what_main_returned(client) -> None:
    r = _post(client)
    assert r.status_code == 200, r.text
    assert _summarise(r.json()) == DEFAULT_PAIRS


def test_the_new_fields_default_to_exactly_tm_conditions(client) -> None:
    """#81's out-of-scope line, asserted rather than trusted.

    Not "some sensible default" — the *same* default, field by field against the
    dict production actually passes to Tm_NN.
    """
    r = _post(client)
    assert r.status_code == 200, r.text
    reported = r.json()[0]["tm_conditions"]
    assert reported["na_mm"] == TM_CONDITIONS["Na"]
    assert reported["k_mm"] == TM_CONDITIONS["K"]
    assert reported["tris_mm"] == TM_CONDITIONS["Tris"]
    assert reported["mg_mm"] == TM_CONDITIONS["Mg"]
    assert reported["dntp_mm"] == TM_CONDITIONS["dNTPs"]
    assert reported["primer_nm"] == TM_CONDITIONS["dnac1"]
    assert reported["template_nm"] == TM_CONDITIONS["dnac2"]


def test_spelling_out_the_defaults_explicitly_changes_nothing(client) -> None:
    """Passing the current values by hand must be indistinguishable from omitting them."""
    r = _post(client, tm_method="nn", na_mm=50, k_mm=0, tris_mm=0, mg_mm=0,
              dntp_mm=0, primer_nm=25, template_nm=25)
    assert r.status_code == 200, r.text
    assert _summarise(r.json()) == DEFAULT_PAIRS


# --- acceptance 3 (reporting): every pair carries its method and conditions ----

def test_every_pair_reports_the_method_and_the_conditions_used(client) -> None:
    r = _post(client)
    assert r.status_code == 200, r.text
    pairs = r.json()
    assert isinstance(pairs, list), "the response must stay a list of pairs"
    assert pairs
    for pair in pairs:
        assert pair["tm_method"] == "nn"
        conditions = pair["tm_conditions"]
        for field in ("na_mm", "k_mm", "tris_mm", "mg_mm", "dntp_mm",
                      "primer_nm", "template_nm"):
            assert field in conditions, field


def test_the_reported_conditions_echo_what_was_asked_for(client) -> None:
    r = _post(client, **PCR_BUFFER)
    assert r.status_code == 200, r.text
    conditions = r.json()[0]["tm_conditions"]
    assert conditions["mg_mm"] == 1.5
    assert conditions["dntp_mm"] == 0.2
    assert conditions["primer_nm"] == 250
    assert conditions["template_nm"] == 0


# --- acceptance 2: tm_method="wallace" actually reaches _tm_wallace -----------

def test_wallace_reports_the_wallace_number_for_every_primer(client) -> None:
    r = _post(client, tm_method="wallace", tm_min=30.0, tm_max=90.0)
    assert r.status_code == 200, r.text
    pairs = r.json()
    assert pairs, "no pairs returned under the Wallace rule"
    for pair in pairs:
        assert pair["tm_method"] == "wallace"
        for role in ("forward", "reverse"):
            p = pair[role]
            assert p["tm"] == pytest.approx(round(_tm_wallace(p["sequence"]), 1)), (
                role, p["sequence"], p["tm"])


def test_wallace_and_nn_actually_disagree_somewhere(client) -> None:
    """Teeth: an option that returned the same numbers would prove nothing.

    Compared per SEQUENCE rather than per result set: the two methods filter on Tm,
    so they legitimately return different primers, and intersecting the two sets
    could quietly come up empty and assert nothing.
    """
    reported = _tm_by_sequence(_post(client, tm_method="wallace",
                                     tm_min=30.0, tm_max=90.0).json())
    assert reported, "no pairs returned under the Wallace rule"
    differing = [seq for seq, tm in reported.items()
                 if abs(tm - _tm_nearest_neighbor(seq)) > 1.0]
    assert differing, (
        "every reported Wallace Tm was within 1 degC of the NN Tm for the same "
        "sequence — this cannot tell the two formulas apart")


# --- acceptance 4: a PCR-like buffer reaches Tm_NN, against recorded values ----

def test_the_pcr_buffer_reproduces_the_recorded_biopython_value(client) -> None:
    """The #79 evidence, at the endpoint: same primer, ~10 degC apart."""
    body = dict(M13_REQUEST)
    default = _tm_by_sequence(client.post("/primers/design", json=body).json())
    body.update(PCR_BUFFER)
    pcr = _tm_by_sequence(client.post("/primers/design", json=body).json())

    assert M13_REVERSE in default, "the M13 primer was not returned under the defaults"
    assert M13_REVERSE in pcr, "the M13 primer was not returned under the PCR buffer"

    assert default[M13_REVERSE] == pytest.approx(M13_TM_DEFAULT_BUFFER, abs=0.05)
    assert pcr[M13_REVERSE] == pytest.approx(M13_TM_PCR_BUFFER, abs=0.05)
    assert pcr[M13_REVERSE] > default[M13_REVERSE]


def test_raising_magnesium_raises_tm(client) -> None:
    """The direction of the Owczarzy correction, as a property over the whole set."""
    default = _tm_by_sequence(_post(client).json())
    salty = _tm_by_sequence(_post(client, mg_mm=1.5, dntp_mm=0.2).json())
    shared = set(default) & set(salty)
    assert shared, "no primer survived both buffers, so this compares nothing"
    assert all(salty[s] > default[s] for s in shared), (
        {s: (default[s], salty[s]) for s in shared if salty[s] <= default[s]})


# --- acceptance 5: unevaluable buffers are refused, not crashed on ------------

def test_no_salt_at_all_is_refused(client) -> None:
    """Tm_NN's salt correction takes log(total ion concentration)."""
    r = _post(client, na_mm=0, k_mm=0, tris_mm=0, mg_mm=0, dntp_mm=0)
    assert 400 <= r.status_code < 500, r.status_code


def test_magnesium_fully_chelated_by_dntps_is_refused(client) -> None:
    """Free Mg2+ is Mg - dNTPs; at or below zero there is no divalent rescue."""
    r = _post(client, na_mm=0, k_mm=0, tris_mm=0, mg_mm=1.5, dntp_mm=1.5)
    assert 400 <= r.status_code < 500, r.status_code


def test_magnesium_alone_is_accepted(client) -> None:
    """The refusal must be about an unevaluable buffer, not about Na being zero."""
    r = _post(client, na_mm=0, k_mm=0, tris_mm=0, mg_mm=1.5, dntp_mm=0.2)
    assert r.status_code == 200, r.text


@pytest.mark.parametrize("primer_nm,template_nm", [(25, 50), (25, 100), (0, 0)])
def test_primer_not_above_half_the_template_is_refused(client, primer_nm, template_nm) -> None:
    """Tm_NN takes log(dnac1 - dnac2/2); at or below zero that is a math domain error."""
    r = _post(client, primer_nm=primer_nm, template_nm=template_nm)
    assert 400 <= r.status_code < 500, r.status_code


@pytest.mark.parametrize("field", ["na_mm", "k_mm", "tris_mm", "mg_mm",
                                   "dntp_mm", "primer_nm", "template_nm"])
def test_negative_concentrations_are_422(client, field) -> None:
    r = _post(client, **{field: -1})
    assert r.status_code == 422, (field, r.status_code)


def test_an_unknown_tm_method_is_422(client) -> None:
    r = _post(client, tm_method="wallis")
    assert r.status_code == 422, r.status_code
