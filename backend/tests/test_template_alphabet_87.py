"""#87 — a template Biopython refuses is refused by both endpoints, with 400.

`Bio.Restriction.FormattedSeq` accepts only `ABCDGHKMNRSTVWY` once layout is
stripped and raises `TypeError("Invalid character found in ...")` on anything
else. So a paste that still has its FASTA `>` header is not an analysable
sequence at all: `/primers/digest` used to let that TypeError out as a 500 while
`/primers/restriction_sites` regex-scanned the raw string and returned 200 with
the header letters counted as bases (EcoRI at 401 instead of 395).

Both are the same missing decision — is this input a sequence? — so both now ask
`backend.sequence_text.validate_template`. This is the #86 sibling and must not
undo it: layout characters still normalise rather than being rejected.
"""

import pytest
from Bio.Restriction import EcoRI
from Bio.Seq import Seq
from fastapi.testclient import TestClient

from backend.sequence_text import normalize_template
from ._puc19 import PUC19

#: pUC19's single EcoRI site, in the real molecule.
ECORI_SITE = 395
#: Where the old regex reported it once ">pUC19" was counted as seven bases.
ECORI_SITE_OVER_HEADER = 401

FASTA_PASTE = ">pUC19\n" + PUC19


def _bodies(template: str) -> list[tuple[str, dict]]:
    body = {"template": template, "enzymes": ["EcoRI"]}
    return [("/primers/digest", body), ("/primers/restriction_sites", body)]


@pytest.mark.parametrize("path", ["/primers/digest", "/primers/restriction_sites"])
def test_fasta_header_is_rejected_with_400_by_both_endpoints(
    client: TestClient, path: str
) -> None:
    """Neither a 500 nor a plausible-looking 200 — a 400 naming the character."""
    r = client.post(path, json={"template": FASTA_PASTE, "enzymes": ["EcoRI"]})
    assert r.status_code == 400, (path, r.status_code, r.text)
    detail = r.json()["detail"]
    assert "'>'" in detail, detail
    assert "position 1" in detail, detail
    # The message is about the user's sequence, not Biopython's internals.
    assert "Invalid character found in" not in detail, detail
    assert str(ECORI_SITE_OVER_HEADER) not in r.text, r.text


def test_the_refusal_matches_the_unknown_enzyme_refusal(client: TestClient) -> None:
    """Same shape as the app's existing bad-input convention, three lines away."""
    unknown = client.post(
        "/primers/digest", json={"template": PUC19, "enzymes": ["NoSuchEnzyme"]}
    )
    assert unknown.status_code == 400
    assert "Unknown enzymes" in unknown.json()["detail"]

    bad_template = client.post(
        "/primers/digest", json={"template": FASTA_PASTE, "enzymes": ["EcoRI"]}
    )
    assert bad_template.status_code == unknown.status_code
    assert isinstance(bad_template.json()["detail"], str)


@pytest.mark.parametrize(
    "junk", [">pUC19\nACGT", "ACGTXACGT", "ACGT-ACGT", "ACGUACGU", "ACGT*", "acgt?"]
)
def test_non_iupac_characters_are_rejected_not_skipped(
    client: TestClient, junk: str
) -> None:
    """A character Biopython refuses is refused — not dropped, not scanned over."""
    for path, body in _bodies(junk):
        r = client.post(path, json=body)
        assert r.status_code == 400, (path, junk, r.status_code, r.text)


def test_layout_characters_are_still_normalised_not_rejected(client: TestClient) -> None:
    """#87's guard must not start rejecting a CRLF or ORIGIN-block paste (#86)."""
    origin = "\n".join(
        f"{i + 1:>9} " + " ".join(PUC19[i + j:i + j + 10].lower() for j in range(0, 60, 10))
        for i in range(0, len(PUC19), 60)
    )
    crlf = "\r\n".join(PUC19[i:i + 60] for i in range(0, len(PUC19), 60))
    for flavour, template in (("origin", origin), ("crlf", crlf), ("tabs", "\tACG\tAATTC\t")):
        for path, body in _bodies(template):
            r = client.post(path, json=body)
            assert r.status_code == 200, (flavour, path, r.status_code, r.text)


def test_iupac_ambiguity_codes_remain_valid_input(client: TestClient) -> None:
    """RE-10's alphabet: N, R, Y and friends are sequence, and keep working."""
    from backend.sequence_text import IUPAC_NUCLEOTIDES

    template = "".join(sorted(IUPAC_NUCLEOTIDES)) + "GAATTC" + "NNNNRYSWKM"
    for path, body in _bodies(template):
        r = client.post(path, json=body)
        assert r.status_code == 200, (path, r.status_code, r.text)
    sites = client.post(
        "/primers/restriction_sites", json={"template": template, "enzymes": ["EcoRI"]}
    ).json()
    assert [s["positions"] for s in sites] == [[len(IUPAC_NUCLEOTIDES)]]


def test_clean_template_is_unaffected(client: TestClient) -> None:
    """The guard costs the good path nothing: pUC19's EcoRI site is still 395."""
    from backend.sequence_text import validate_template

    sites = client.post(
        "/primers/restriction_sites", json={"template": PUC19, "enzymes": ["EcoRI"]}
    ).json()
    assert [s["positions"] for s in sites] == [[ECORI_SITE]]
    assert validate_template(PUC19) == PUC19


def test_validate_template_agrees_with_biopython_on_every_ascii_character() -> None:
    """The alphabet is Bio.Restriction's, checked against it rather than assumed.

    IUPACData.ambiguous_dna_values lists a legacy "X", which FormattedSeq's table
    does not accept — validating against that set would pass an input that then
    crashes, so this pins agreement with the code that does the searching.
    """
    from fastapi import HTTPException

    from backend.sequence_text import validate_template

    for code in range(32, 127):
        char = chr(code)
        probe = "ACGT" + char + "ACGT"
        try:
            EcoRI.search(Seq(probe), linear=True)
        except TypeError:
            refused_by_biopython = True
        else:
            refused_by_biopython = False

        if normalize_template(probe) == "ACGTACGT":
            continue  # layout: stripped by both, nothing to disagree about
        try:
            validate_template(probe)
        except HTTPException as e:
            assert refused_by_biopython, (char, e.detail)
        else:
            assert not refused_by_biopython, char


def test_first_non_iupac_reports_the_offending_position() -> None:
    """The index is into the molecule as the app reads it, layout already gone."""
    from backend.sequence_text import first_non_iupac

    assert first_non_iupac("ACGTN") is None
    assert first_non_iupac("") is None
    assert first_non_iupac(">PUCACGT") == (0, ">")
    assert first_non_iupac(normalize_template("acgt acgt x")) == (8, "X")
