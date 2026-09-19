"""#86 — a pasted template's layout characters must not move its coordinates.

A GenBank ORIGIN block and a CRLF-wrapped FASTA body denote the same molecule as
the bare sequence. `Bio.Restriction.FormattedSeq` strips `string.whitespace +
string.digits` before searching, and `/primers/digest` reports the positions that
search returns — so anything the app strips *less* than that leaves it slicing a
string Biopython never saw, while `/primers/restriction_sites` indexes the dirty
string with its own regex. The two endpoints then disagree about the same
request body.

The sealed gx suite pins this against the real pUC19; these are the same
properties on the app's own fixtures, so `make check` catches a regression.
"""

from fastapi.testclient import TestClient

from backend.sequence_text import normalize_template
from ._puc19 import PUC19

#: pUC19's single EcoRI site and the cut one base into it (EcoRI fst5 == 1).
ECORI_SITE = 395
ECORI_CUT = 396


def genbank_origin(seq: str) -> str:
    """The sequence as it comes out of a GenBank ORIGIN block."""
    lines = []
    for i in range(0, len(seq), 60):
        blocks = " ".join(seq[i + j:i + j + 10].lower() for j in range(0, 60, 10))
        lines.append(f"{i + 1:>9} {blocks}")
    return "\n".join(lines)


def crlf(seq: str) -> str:
    return "\r\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))


def test_normalize_template_strips_layout_but_not_letters() -> None:
    """Whitespace and digits go; anything else — including junk — stays."""
    assert normalize_template("  ac g\r\nt\t") == "ACGT"
    assert normalize_template("   1 acgtac gtac\n  61 gtacgtacgt") == "ACGTACGTACGTACGTACGT"
    # Not this issue's job to reject a FASTA header — that is #87. The digits in
    # one go the way any other digits do; what survives is the non-IUPAC junk
    # that makes it a rejection case rather than a normalisation one.
    assert normalize_template(">pUC19\nACGT") == ">PUCACGT"


def test_restriction_sites_ignore_line_numbers_and_carriage_returns(client: TestClient) -> None:
    """RE-11: the reported site is the site in the molecule, whatever the layout."""
    for flavour, template in (("genbank_origin", genbank_origin(PUC19)), ("crlf", crlf(PUC19))):
        r = client.post("/primers/restriction_sites",
                        json={"template": template, "enzymes": ["EcoRI"]})
        assert r.status_code == 200, flavour
        assert [s["positions"] for s in r.json()] == [[ECORI_SITE]], flavour


def test_digest_reports_the_molecule_not_the_paste(client: TestClient) -> None:
    """RE-11: template_length, cuts and fragment lengths are the clean molecule's."""
    for flavour, template in (("genbank_origin", genbank_origin(PUC19)), ("crlf", crlf(PUC19))):
        r = client.post("/primers/digest", json={"template": template, "enzymes": ["EcoRI"]})
        assert r.status_code == 200, flavour
        result = r.json()
        assert result["template_length"] == len(PUC19) == 2686, flavour
        assert result["cut_positions"] == [ECORI_CUT], flavour
        assert sorted(f["length"] for f in result["fragments"]) == [396, 2290], flavour


def test_digest_fragment_sequences_are_nucleotides(client: TestClient) -> None:
    """RE-11: `DigestFragment.sequence` is DNA — the Assembly tab ligates it."""
    for flavour, template in (("genbank_origin", genbank_origin(PUC19)), ("crlf", crlf(PUC19))):
        r = client.post("/primers/digest", json={"template": template, "enzymes": ["EcoRI"]})
        for frag in r.json()["fragments"]:
            stray = sorted({c for c in frag["sequence"] if c not in "ACGT"})
            assert stray == [], (flavour, stray)
            assert len(frag["sequence"]) == frag["length"], flavour


def test_both_endpoints_share_one_coordinate_system(client: TestClient) -> None:
    """RE-12: site + fst5 == cut, for the same request body, EcoRI fst5 == 1."""
    for flavour, template in (("genbank_origin", genbank_origin(PUC19)),
                              ("crlf", crlf(PUC19)),
                              ("clean", PUC19)):
        body = {"template": template, "enzymes": ["EcoRI"]}
        starts = [p for s in client.post("/primers/restriction_sites", json=body).json()
                  for p in s["positions"]]
        cuts = client.post("/primers/digest", json=body).json()["cut_positions"]
        assert [p + 1 for p in starts] == cuts, flavour


def test_site_straddling_a_line_break_is_not_silently_missed(client: TestClient) -> None:
    """A 60-column wrap puts a break inside a site sooner or later (#86)."""
    probe = "AAG\r\nAATTCAA"
    r = client.post("/primers/restriction_sites", json={"template": probe, "enzymes": ["EcoRI"]})
    assert [s["positions"] for s in r.json()] == [[2]]
    d = client.post("/primers/digest", json={"template": probe, "enzymes": ["EcoRI"]})
    assert d.json()["cut_positions"] == [3]
    assert d.json()["template_length"] == 10


def test_design_primers_normalises_the_template_too(client: TestClient) -> None:
    """The same box feeds PCR design; its product sizes are the molecule's."""
    region = PUC19[:400]
    body = {"product_min": 100, "product_max": 300, "primer_len_min": 18,
            "primer_len_max": 22, "tm_min": 50.0, "tm_max": 70.0, "max_pairs": 3}
    clean = client.post("/primers/design", json={"template": region, **body}).json()
    dirty = client.post("/primers/design", json={"template": crlf(region), **body}).json()
    assert dirty == clean
