"""Shared fixtures and INDEPENDENT oracle helpers for the gx restriction suite.

Nothing in here imports the subject's own search/pattern helpers. The IUPAC
expansion used by the oracles is built from the pinned corpus copy of
Biopython's ``Bio.Data.IUPACData.ambiguous_dna_values``
(gx/corpus/biopython-1.85/Bio/Data/IUPACData.py), which is the same table the
REBASE-derived ``compsite`` regexes in Restriction_Dictionary.py encode.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest
from Bio.Data.IUPACData import ambiguous_dna_values

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

CORPUS = Path(__file__).resolve().parent / "corpus"
PUC19_PATH = CORPUS / "inputs" / "puc19_L09137.2.txt"

_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
# Deterministic preference when several codes denote the same base set (X/N).
_CODE_ORDER = "ACGTRYSWKMBDHVN"


def iupac_regex(site: str) -> str:
    """Regex for a recognition sequence, expanding ambiguity codes by hand."""
    out = []
    for ch in site.upper():
        bases = ambiguous_dna_values[ch]
        out.append(bases if len(bases) == 1 else "[" + "".join(sorted(bases)) + "]")
    return "".join(out)


def iupac_complement(code: str) -> str:
    """The ambiguity code denoting the complement of ``code``'s base set."""
    want = {_COMPLEMENT[b] for b in ambiguous_dna_values[code.upper()]}
    for candidate in _CODE_ORDER:
        if set(ambiguous_dna_values[candidate]) == want:
            return candidate
    raise AssertionError(f"no IUPAC code complements {code!r}")


def iupac_reverse_complement(site: str) -> str:
    return "".join(iupac_complement(c) for c in reversed(site.upper()))


def find_starts(template: str, site: str) -> list[int]:
    """0-based starts of every (possibly overlapping) occurrence of ``site``."""
    return [m.start() for m in re.finditer(f"(?={iupac_regex(site)})", template)]


def find_starts_circular(template: str, site: str) -> list[int]:
    """0-based starts on a circle: every window of the doubled sequence."""
    n, k = len(template), len(site)
    if n < k:
        return []
    doubled = template + template
    return [i for i in range(n) if re.match(iupac_regex(site), doubled[i : i + k])]


@pytest.fixture(scope="session")
def puc19() -> str:
    """The pinned 2686 bp pUC19 test input (gx/corpus/inputs, see MANIFEST.json)."""
    seq = PUC19_PATH.read_text().strip()
    assert len(seq) == 2686, len(seq)
    assert set(seq) <= set("ACGT")
    return seq


@pytest.fixture
def panel_with(monkeypatch):
    """Add an enzyme to the panel dicts the endpoints read, then restore.

    Used to drive the subject's own site finder with a recognition sequence that
    is not its own reverse complement (RE-3): the shipped panel happens to be
    entirely palindromic (RE-2), so a top-strand-only scan cannot be exercised
    through the shipped enzyme list alone.
    """
    from Bio import Restriction

    from backend.routers import primers

    def _add(name: str) -> None:
        enz = getattr(Restriction, name)
        monkeypatch.setitem(primers.RESTRICTION_ENZYMES, name, enz.site)
        monkeypatch.setitem(primers._ENZYMES, name, enz)

    return _add
