"""Fixtures for the gx CRISPR round. Oracles live in ``crispr_oracle.py``.

Nothing is imported from here by name -- see the note at the top of
``crispr_oracle``: ``from conftest import ...`` is what forced earlier rounds
into separate processes, so this file is fixtures only.

The subject is ``backend/routers/crispr.py`` (``design_guides`` /
``score_guide``) together with ``backend/crispr_offtarget.py``, called
in-process exactly as ``gui/panels/crispr_panel.py:_run_design`` calls them.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from Bio import SeqIO

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

INPUTS = REPO_ROOT / "gx" / "corpus" / "inputs"


def _genbank_sequence(path: Path) -> str:
    return str(SeqIO.read(path, "genbank").seq).upper()


@pytest.fixture(scope="session")
def puc19() -> str:
    """pUC19, 2686 bp, circular plasmid (gx/corpus/inputs, see MANIFEST.json)."""
    seq = INPUTS / "puc19_L09137.2.txt"
    text = seq.read_text().strip()
    assert len(text) == 2686 and set(text) <= set("ACGT")
    return text


@pytest.fixture(scope="session")
def jcv() -> str:
    """NC_001699.1, JC polyomavirus, 5130 bp, circular genome."""
    seq = _genbank_sequence(INPUTS / "NC_001699.1.gb")
    assert len(seq) == 5130
    return seq


@pytest.fixture(scope="session")
def bfdv() -> str:
    """AF071878.1, beak and feather disease virus, 1993 bp, circular genome."""
    seq = _genbank_sequence(INPUTS / "AF071878.1.gb")
    assert len(seq) == 1993
    return seq


@pytest.fixture(scope="session")
def col1a1() -> str:
    """NG_007400.1, human COL1A1 RefSeqGene, 24,544 bp -- the real-scale input."""
    seq = _genbank_sequence(INPUTS / "NG_007400.1.gb")
    assert len(seq) == 24544 and set(seq) <= set("ACGT")
    return seq


@pytest.fixture(scope="session")
def corpus(puc19, jcv, bfdv, col1a1) -> dict[str, str]:
    """Every pinned full-length input, by the name the register uses."""
    return {
        "puc19_L09137.2": puc19,
        "NC_001699.1": jcv,
        "AF071878.1": bfdv,
        "NG_007400.1": col1a1,
    }


@pytest.fixture
def design():
    """Call the subject's design endpoint the way the GUI worker does."""
    from backend.models.crispr import CasType, CRISPRDesignRequest
    from backend.routers.crispr import design_guides

    def _design(target: str, cas: str, **kw):
        return design_guides(
            CRISPRDesignRequest(target_sequence=target, cas_type=CasType(cas), **kw)
        )

    return _design
