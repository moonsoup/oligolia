"""Differential oracle for protein properties — #63, second piece of #75's layer.

Four independent defects were confirmed against Biopython's ProtParam:

  extinction coefficient  cystine term double-counted: `(C//2) * 125 * 2`
  instability index       25 of 400 DIWV entries shipped; every missing pair
                          silently contributed 0.0, so proteins read as stable
  molecular weight        monoisotopic residue masses reported as average Da
  isoelectric point       a pKa set described as Bjellqvist that is not

The oracle is `Bio.SeqUtils.ProtParam.ProteinAnalysis`, already a dependency and
carrying the real 400-entry DIWV table (`Bio.SeqUtils.ProtParamData.DIWV`).

WHAT IS AND IS NOT INDEPENDENT HERE. Production now *wraps* ProtParam, so a live
ProtParam call in the test would compare a function to itself — the same hit Codex
made against the Tm file. So:

  * `PROTEIN_REFERENCE` holds values RECORDED from Biopython 1.85 on 2026-09-17.
    Nothing computes them at run time, so a change in the wrapper, in Biopython's
    tables, or in the substitution map all surface as a failure. This is the pin.
  * live-ProtParam comparisons are kept and labelled as wrapper-integrity checks.
  * the independent assertions are the ones comparing DIFFERENT formulas: the old
    doubled cystine term against the correct one, and cross-protein distinctness
    of the instability index.

The recorded values independently corroborate the audit in the issue: lysozyme
instability 16.09 and pI 9.32, HBB oxidised extinction 15595, lysozyme 37970 —
all reproduced here from a separate sequence transcription.
"""

from __future__ import annotations

import pytest
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Real sequences, no truncation. HBB is the canonical human beta-globin chain.
HBB = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVL"
    "GAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGV"
    "ANALAHKYH"
)
# Hen egg white lysozyme, the other protein the audit used.
LYSOZYME = (
    "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRT"
    "PGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
)
# A cysteine-free protein, so the extinction path is exercised both ways.
NO_CYS = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAA"

PROTEINS = {"HBB": HBB, "lysozyme": LYSOZYME, "no_cys": NO_CYS}


#: Recorded from Biopython 1.85 on 2026-09-17. Regenerate deliberately.
PROTEIN_REFERENCE = {
    "HBB": dict(mw=15998.2064, pi=6.7445, ii=6.1578, ext_red=15470, ext_ox=15595, gravy=0.0136),
    "lysozyme": dict(mw=14313.004, pi=9.3238, ii=16.0915, ext_red=37470, ext_ox=37970, gravy=-0.4721),
    "no_cys": dict(mw=7852.5465, pi=4.05, ii=92.8386, ext_red=11000, ext_ox=11000, gravy=-0.7657),
}


def oracle(seq: str) -> ProteinAnalysis:
    """A live ProtParam object. Wrapper-integrity only — see the module docstring."""
    return ProteinAnalysis(seq)


@pytest.mark.parametrize("name", sorted(PROTEIN_REFERENCE))
def test_all_properties_match_the_recorded_reference(client, name: str) -> None:
    """The real pin: recorded values, not a live call to the wrapped library."""
    ref = PROTEIN_REFERENCE[name]
    got = _props(client, PROTEINS[name])
    assert abs(got["molecular_weight_da"] - ref["mw"]) <= 1.0, (name, got["molecular_weight_da"], ref["mw"])
    assert abs(got["isoelectric_point"] - ref["pi"]) <= 0.1, (name, got["isoelectric_point"], ref["pi"])
    assert abs(got["instability_index"] - ref["ii"]) <= 0.1, (name, got["instability_index"], ref["ii"])
    assert abs(got["extinction_coeff_red"] - ref["ext_red"]) <= 1, (name, got["extinction_coeff_red"], ref["ext_red"])
    assert abs(got["extinction_coeff_ox"] - ref["ext_ox"]) <= 1, (name, got["extinction_coeff_ox"], ref["ext_ox"])
    assert abs(got["gravy"] - ref["gravy"]) <= 0.01, (name, got["gravy"], ref["gravy"])


def test_the_recorded_reference_still_agrees_with_this_biopython() -> None:
    """If ProtParam's tables move, fail loudly rather than following silently."""
    drifted = {}
    for name, seq in PROTEINS.items():
        a = oracle(seq)
        red, ox = a.molar_extinction_coefficient()
        live = dict(mw=a.molecular_weight(), pi=a.isoelectric_point(),
                    ii=a.instability_index(), ext_red=int(red), ext_ox=int(ox),
                    gravy=a.gravy())
        ref = PROTEIN_REFERENCE[name]
        for k, v in ref.items():
            if abs(live[k] - v) > 0.01:
                drifted[f"{name}.{k}"] = (v, round(live[k], 4))
    assert not drifted, (
        "Biopython's ProtParam no longer matches the values recorded on 2026-09-17; "
        f"re-record deliberately and note why: {drifted}"
    )


def _props(client, seq: str) -> dict:
    r = client.post("/analysis/protein_properties", params={"sequence": seq})
    assert r.status_code == 200, r.text
    return r.json()


@pytest.mark.parametrize("name", sorted(PROTEINS))
def test_molecular_weight_is_average_mass(client, name: str) -> None:
    """Wrapper integrity. Reported as Da, so it must be average, not monoisotopic."""
    seq = PROTEINS[name]
    got = _props(client, seq)["molecular_weight_da"]
    expected = oracle(seq).molecular_weight()
    # Average vs monoisotopic differ by ~0.06% — well outside 1 Da on a 15 kDa chain.
    assert abs(got - expected) <= 1.0, (name, got, expected)


@pytest.mark.parametrize("name", sorted(PROTEINS))
def test_instability_index_matches_protparam(client, name: str) -> None:
    """Wrapper integrity — not independent; production wraps this function."""
    seq = PROTEINS[name]
    got = _props(client, seq)["instability_index"]
    expected = oracle(seq).instability_index()
    assert abs(got - expected) <= 0.1, (name, got, expected)


def test_the_instability_index_is_not_near_zero_for_everything() -> None:
    """The shape of the old bug: 25 of 400 dipeptides scored, the rest 0.0.

    A table that covers 6% of dipeptides makes every protein look stable, which is
    the clinically misleading direction.
    """
    values = [oracle(s).instability_index() for s in PROTEINS.values()]
    assert max(values) > 5.0, values
    assert len({round(v, 1) for v in values}) == len(values), (
        "different proteins should not share an instability index"
    )


@pytest.mark.parametrize("name", sorted(PROTEINS))
def test_extinction_coefficients_match_protparam(client, name: str) -> None:
    """Gill & von Hippel: 125 is per CYSTINE, so (C//2)*125 — not times two."""
    seq = PROTEINS[name]
    props = _props(client, seq)
    reduced, oxidised = oracle(seq).molar_extinction_coefficient()
    assert props["extinction_coeff_red"] == pytest.approx(reduced, abs=1), (name, props)
    assert props["extinction_coeff_ox"] == pytest.approx(oxidised, abs=1), (name, props)


def test_the_old_cystine_term_really_was_double() -> None:
    """Teeth: the old formula must differ from the oracle on a cysteine protein."""
    seq = LYSOZYME
    c = seq.count("C")
    assert c >= 2, "this test needs a protein with cystines"
    old_ox = seq.count("W") * 5500 + seq.count("Y") * 1490 + (c // 2) * 125 * 2
    _reduced, oxidised = oracle(seq).molar_extinction_coefficient()
    assert old_ox > oxidised, (old_ox, oxidised)
    assert old_ox - oxidised == pytest.approx((c // 2) * 125, abs=1), (
        "the excess should be exactly one extra cystine term"
    )


@pytest.mark.parametrize("name", sorted(PROTEINS))
def test_isoelectric_point_matches_protparam(client, name: str) -> None:
    seq = PROTEINS[name]
    got = _props(client, seq)["isoelectric_point"]
    expected = oracle(seq).isoelectric_point()
    assert abs(got - expected) <= 0.1, (name, got, expected)


@pytest.mark.parametrize("name", sorted(PROTEINS))
def test_gravy_matches_protparam(client, name: str) -> None:
    """Kyte-Doolittle GRAVY — the one the old code already got right; pin it."""
    seq = PROTEINS[name]
    got = _props(client, seq)["gravy"]
    expected = oracle(seq).gravy()
    assert abs(got - expected) <= 0.01, (name, got, expected)


# --- non-standard residues must stay supported, and the substitution declared ---

def test_selenocysteine_is_accepted_and_the_substitution_is_declared(client) -> None:
    """Delegating to ProtParam broke this: it raises KeyError on U.

    The endpoint accepts selenocysteine and pyrrolysine, so refusing them would
    have removed a working capability. Substituting the chemical analogue is fine;
    doing it silently is not, which is why it comes back in the response.
    """
    r = client.post("/analysis/protein_properties", params={"sequence": "MACUGPWK"})
    assert r.status_code == 200, r.text
    data = r.json()
    # The real residue is still reported.
    assert data["aa_composition"].get("U", 0) == 1
    # And the substitution is on the record.
    assert data["nonstandard_substituted"] == {"U": "C"}
    assert data["molecular_weight_da"] > 0
    assert data["instability_index"] != 0.0


def test_an_all_standard_sequence_declares_no_substitution(client) -> None:
    r = client.post("/analysis/protein_properties", params={"sequence": HBB})
    assert r.status_code == 200, r.text
    assert r.json()["nonstandard_substituted"] == {}


@pytest.mark.parametrize("code,analogue", [("U", "C"), ("O", "K"), ("B", "D"), ("Z", "E"), ("J", "L"), ("X", "G")])
def test_every_nonstandard_code_is_handled_and_named(client, code: str, analogue: str) -> None:
    r = client.post("/analysis/protein_properties", params={"sequence": f"MACG{code}PWK"})
    assert r.status_code == 200, (code, r.text)
    assert r.json()["nonstandard_substituted"] == {code: analogue}
