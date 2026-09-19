"""INDEPENDENT oracles for the coordinates round-2 suite (CO2-1 … CO2-3).

Round 2 corrects the sealed round-1 test
``gx/coords/test_coords_edit_endpoint.py::test_co_14_surviving_annotations_point_at_the_right_bases``
(``gx-co-53``); see ``SUPERSEDES.md`` for the defect and the evidence.

Deliberately NOT a ``conftest.py``: round 1 already ships
``gx/conftest.py`` and ``gx/coords/conftest.py``, and a round-2 test that said
``from conftest import ...`` would bind to whichever of those landed in
``sys.path`` first. Import this module by its own name instead::

    from coords_r2_oracle import requirement, expected_extract

Nothing here calls the subject. Every expected value comes from one of three
places, and never from the app:

* ``Bio.SeqIO`` / ``SeqFeature.extract`` run on the pinned NC_001699.1 bytes —
  the location semantics pinned in ``gx/corpus/biopython-1.85/Bio/SeqFeature.py``;
* a re-derivation by index from the raw sequence string (``seq[s:e]``);
* a literal recorded in ``gx/coords/r2/register.json``.
"""

from __future__ import annotations

import hashlib
import json
import sys
from io import StringIO
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

JCV = REPO_ROOT / "gx" / "corpus" / "inputs" / "NC_001699.1.gb"
REGISTER = json.loads((HERE / "register.json").read_text())


def requirement(rid: str) -> dict:
    """The pinned requirement, so tests read their numbers from the register."""
    for r in REGISTER:
        if r["id"] == rid:
            return r
    raise AssertionError(f"no requirement {rid!r} in gx/coords/r2/register.json")


# ── oracle: Biopython on the pinned bytes ────────────────────────────────────

def jcv_text() -> str:
    return JCV.read_text()


def jcv_bio() -> SeqRecord:
    rec = SeqIO.read(StringIO(jcv_text()), "genbank")
    assert len(rec.seq) == 5130 and len(rec.features) == 28
    return rec


def bio_parts(feature) -> list[tuple[int, int]]:
    """Every interval the feature occupies, in Biopython's 5'-to-3' part order."""
    return [(int(p.start), int(p.end)) for p in feature.location.parts]


def bio_strand(feature) -> str:
    return {1: "+", -1: "-"}.get(feature.location.strand, ".")


def bio_extract(feature, record: SeqRecord) -> str:
    return str(feature.extract(record.seq))


# ── oracle: re-derivation by index ───────────────────────────────────────────

def spliced(seq: str, parts: list[tuple[int, int]], strand: str) -> str:
    """The subsequence a feature denotes, built straight out of the string.

    Parts are concatenated in the order given (INSDC join semantics); on the
    minus strand each part is reverse-complemented, matching
    ``SeqFeature.extract``.
    """
    pieces = [seq[s:e] for s, e in parts]
    if strand == "-":
        return "".join(str(Seq(p).reverse_complement()) for p in pieces)
    return "".join(pieces)


def contains(part: tuple[int, int], start: int, end: int) -> bool:
    """Does this part strictly contain the edit ``[start, end)``?

    Strictly: a part flush with either boundary does not contain the edit — the
    inserted bases fall outside it, and a deletion that reaches its edge takes
    bases it does not own.
    """
    p_start, p_end = part
    if start == end:
        return p_start < start < p_end
    return p_start <= start and p_end >= end and not (p_end <= start or p_start >= end)


def overlaps_partially(part: tuple[int, int], start: int, end: int) -> bool:
    """Does the edit take some but not all of this part's claim?"""
    p_start, p_end = part
    if start == end:
        return False
    return p_start < end and p_end > start and not contains(part, start, end)


def expected_parts(
    parts: list[tuple[int, int]], start: int, end: int, delta: int
) -> list[tuple[int, int]]:
    """Where each part must land, derived from the edit alone.

    A part strictly before the edit does not move; one strictly after it moves
    by ``delta``; one that contains it keeps its start and moves its end.
    """
    out = []
    for p_start, p_end in parts:
        if p_end <= start:
            out.append((p_start, p_end))
        elif p_start >= end:
            out.append((p_start + delta, p_end + delta))
        elif contains((p_start, p_end), start, end):
            out.append((p_start, p_end + delta))
        else:
            raise AssertionError(
                f"part ({p_start}, {p_end}) partially overlaps [{start}, {end}); "
                "no round-2 requirement pins that case"
            )
    return out


def expected_extract(
    parts: list[tuple[int, int]],
    strand: str,
    original_seq: str,
    start: int,
    end: int,
    inserted: str,
) -> str:
    """The bases the feature must read after the edit, built from the ORIGINAL string.

    Carried parts contribute their own bases verbatim. A part that contains the
    edit contributes its bases with ``[start, end)`` swapped for ``inserted`` —
    those inserted bases are now inside its claim. Nothing here consults the
    edited record the subject produced.
    """
    pieces = []
    for p_start, p_end in parts:
        if p_end <= start or p_start >= end:
            pieces.append(original_seq[p_start:p_end])
        elif contains((p_start, p_end), start, end):
            pieces.append(
                original_seq[p_start:start] + inserted + original_seq[end:p_end]
            )
        else:
            raise AssertionError(
                f"part ({p_start}, {p_end}) partially overlaps [{start}, {end})"
            )
    if strand == "-":
        return "".join(str(Seq(p).reverse_complement()) for p in pieces)
    return "".join(pieces)


def sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()
