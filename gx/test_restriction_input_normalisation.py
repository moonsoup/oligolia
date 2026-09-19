"""Template normalisation and cross-endpoint coherence: RE-11, RE-12, RE-13.

Biopython's ``FormattedSeq`` — which is what ``/primers/digest`` delegates its
cut positions to — strips ``string.whitespace + string.digits`` from the
sequence before searching, and refuses anything left that is not an IUPAC
letter. So a template pasted with layout characters (a GenBank ORIGIN block, or
CRLF line endings) denotes the same molecule as the clean sequence, and both
endpoints must describe that molecule in one coordinate system.

Oracle: the clean 2686 bp pUC19 and ``Bio.Restriction`` on it.
"""

from __future__ import annotations

import pytest
from Bio import Restriction
from Bio.Seq import Seq
from fastapi.testclient import TestClient

from backend.main import app
from backend.routers.primers import (
    DigestRequest,
    RestrictionRequest,
    digest,
    restriction_sites,
)


def genbank_origin_block(seq: str) -> str:
    """pUC19 as pasted out of a GenBank record: 60 bases a line, numbered."""
    lines = []
    for i in range(0, len(seq), 60):
        blocks = " ".join(seq[i + j : i + j + 10].lower() for j in range(0, 60, 10))
        lines.append(f"{i + 1:>9} {blocks}")
    return "\n".join(lines)


def crlf_wrapped(seq: str) -> str:
    return "\r\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))


@pytest.fixture
def templates(puc19):
    return {
        "clean": puc19,
        "genbank_origin": genbank_origin_block(puc19),
        "crlf": crlf_wrapped(puc19),
    }


@pytest.mark.parametrize("flavour", ["clean", "genbank_origin", "crlf"])
def test_re_11_layout_characters_do_not_move_the_coordinates(templates, puc19, flavour):
    """RE-11: site 395, cut 396, length 2686, fragments 396 + 2290, whatever the layout."""
    template = templates[flavour]
    # The oracle is the clean molecule: one EcoRI site at 395, cut at 396.
    assert Restriction.EcoRI.search(Seq(puc19), linear=True) == [397]

    sites = restriction_sites(RestrictionRequest(template=template, enzymes=["EcoRI"]))
    assert [s.positions for s in sites] == [[395]]

    result = digest(DigestRequest(template=template, enzymes=["EcoRI"]))
    assert result.template_length == 2686
    assert result.cut_positions == [396]
    assert sorted(f.length for f in result.fragments) == [396, 2290]
    assert sum(f.length for f in result.fragments) == 2686


@pytest.mark.parametrize("flavour", ["clean", "genbank_origin", "crlf"])
def test_re_11_fragment_sequences_are_nucleotides_only(templates, flavour):
    """RE-11: a reported fragment sequence is DNA — no digits, no control characters."""
    result = digest(DigestRequest(template=templates[flavour], enzymes=["EcoRI"]))
    for frag in result.fragments:
        stray = sorted({c for c in frag.sequence if c not in "ACGT"})
        assert stray == [], stray


@pytest.mark.parametrize("flavour", ["clean", "genbank_origin", "crlf"])
def test_re_12_site_and_cut_endpoints_share_one_coordinate_system(templates, flavour):
    """RE-12: len(sites) == len(cuts) and cut == site start + fst5, same request body."""
    template = templates[flavour]
    fst5 = Restriction.EcoRI.fst5
    assert fst5 == 1

    sites = restriction_sites(RestrictionRequest(template=template, enzymes=["EcoRI"]))
    starts = [p for s in sites for p in s.positions]
    cuts = digest(DigestRequest(template=template, enzymes=["EcoRI"])).cut_positions

    assert len(starts) == len(cuts)
    assert [p + fst5 for p in starts] == cuts


def test_re_12_coherence_holds_across_the_panel(puc19):
    """RE-12: the coherence rule holds for every panel enzyme on the clean sequence."""
    from backend.routers.primers import _CURATED_ENZYMES

    for name in _CURATED_ENZYMES:
        fst5 = getattr(Restriction, name).fst5
        starts = [
            p
            for s in restriction_sites(RestrictionRequest(template=puc19, enzymes=[name]))
            for p in s.positions
        ]
        cuts = digest(DigestRequest(template=puc19, enzymes=[name])).cut_positions
        assert [p + fst5 for p in starts] == cuts, name


def test_re_13_non_iupac_template_is_rejected_by_both_endpoints(puc19):
    """RE-13: a FASTA header left in the paste is not a sequence — reject it, 4xx."""
    template = ">pUC19\n" + puc19
    with pytest.raises(TypeError, match="Invalid character found"):
        # Biopython refuses the raw string, so it is not an analysable sequence.
        Restriction.EcoRI.search(Seq(template), linear=True)

    client = TestClient(app, raise_server_exceptions=False)

    digest_response = client.post(
        "/primers/digest", json={"template": template, "enzymes": ["EcoRI"]}
    )
    assert 400 <= digest_response.status_code < 500, digest_response.status_code

    sites_response = client.post(
        "/primers/restriction_sites", json={"template": template, "enzymes": ["EcoRI"]}
    )
    assert 400 <= sites_response.status_code < 500, (
        sites_response.status_code,
        sites_response.json(),
    )


def test_re_13_unknown_enzyme_is_the_reference_rejection(puc19):
    """RE-13: the app's own bad-input convention is HTTP 400 with a detail message."""
    client = TestClient(app, raise_server_exceptions=False)
    response = client.post(
        "/primers/digest", json={"template": puc19, "enzymes": ["NoSuchEnzyme"]}
    )
    assert response.status_code == 400
    assert "Unknown enzymes" in response.json()["detail"]
