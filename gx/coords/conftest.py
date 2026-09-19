"""Fixtures and INDEPENDENT oracles for the gx coordinate-fidelity suite.

The subject is the app's record pipeline: ``backend/formats/genbank.py``
(``read_genbank`` / ``read_embl`` / ``write_genbank``),
``backend/services/annotations.py`` (``shift_annotations`` / ``flip_annotations``,
the arithmetic ``gui/panels/sequence_panel.py:_commit_edit`` calls on every edit)
and ``backend/routers/sequences.py``'s edit endpoint.

Nothing here calls any of those. Every expected value comes from one of three
places, and never from the app:

* ``Bio.SeqIO`` / ``SeqFeature.extract`` run on the same pinned bytes — the
  location semantics pinned in ``gx/corpus/biopython-1.85/Bio/SeqFeature.py``,
  ``Bio/SeqIO/InsdcIO.py`` and ``Bio/GenBank/__init__.py``;
* a re-derivation by index from the raw sequence string (``seq[s:e]``);
* a literal recorded in ``gx/coords/register.json``.

The INSDC Feature Table Definition (``gx/corpus/insdc/FT_current.txt``) is the
authority for what a location *means* — join vs order, ``<``/``>``, remote
references, and origin-spanning joins on circular records.
"""

from __future__ import annotations

import json
import sys
from io import StringIO
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

CORPUS = REPO_ROOT / "gx" / "corpus"
INPUTS = CORPUS / "inputs"
REGISTER = json.loads((Path(__file__).resolve().parent / "register.json").read_text())

JCV = INPUTS / "NC_001699.1.gb"          # 5130 bp, circular, origin-spanning join,
                                          # complement(join(...)), complement(<...)
BFDV_GB = INPUTS / "AF071878.1.gb"       # 1993 bp, circular; NCBI's copy
BFDV_EMBL = INPUTS / "AF071878.1.embl"   # the same INSDC record, EBI's copy


def requirement(rid: str) -> dict:
    """The pinned requirement, so tests read their numbers from the register."""
    for r in REGISTER:
        if r["id"] == rid:
            return r
    raise AssertionError(f"no requirement {rid!r} in gx/coords/register.json")


# ── oracle: Biopython on the same bytes ──────────────────────────────────────

def bio_record(path: Path, fmt: str = "genbank") -> SeqRecord:
    return SeqIO.read(StringIO(path.read_text()), fmt)


def bio_parts(feature) -> list[tuple[int, int]]:
    """Every interval the feature occupies, in Biopython's 5'-to-3' part order."""
    return [(int(p.start), int(p.end)) for p in feature.location.parts]


def bio_strand(feature) -> str:
    return {1: "+", -1: "-"}.get(feature.location.strand, ".")


def bio_extract(feature, record: SeqRecord) -> str:
    return str(feature.extract(record.seq))


def bio_profile(record: SeqRecord) -> list[tuple]:
    """(type, parts, strand, extract) per feature — the whole oracle in one shot."""
    return [
        (f.type, bio_parts(f), bio_strand(f), bio_extract(f, record))
        for f in record.features
    ]


def reread(genbank_text: str) -> SeqRecord:
    """Parse the subject's own export with the authority, not with the subject."""
    return SeqIO.read(StringIO(genbank_text), "genbank")


def feature_lines(genbank_text: str) -> list[tuple[str, str]]:
    """(key, location) for every top-level line of an exported FEATURES block."""
    out, in_features = [], False
    for line in genbank_text.splitlines():
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if in_features:
            if line.startswith(("ORIGIN", "//")) or (line[:1] not in (" ", "")):
                break
            if line.startswith(" " * 5) and line[5:6] != " ":
                key, _, loc = line[5:].strip().partition(" ")
                out.append((key, loc.strip()))
    return out


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


def touched(parts: list[tuple[int, int]], start: int, end: int) -> bool:
    """Does an edit replacing [start, end) change any of these bases?

    An insertion (start == end) only touches a part if it lands strictly inside
    it; at a boundary the inserted bases are outside the feature either way.
    """
    if start == end:
        return any(p_start < start < p_end for p_start, p_end in parts)
    return any(p_start < end and p_end > start for p_start, p_end in parts)


# ── building an edited record through the subject's own arithmetic ───────────

def edited_sequence(seq, *, start: int, end: int, inserted: str):
    """Apply an edit the way ``sequence_panel._commit_edit`` does.

    Returns ``(new_sequence_model, dropped_annotations)``. The splice and the
    ``shift_annotations`` call are exactly the pair the GUI makes at
    ``gui/panels/sequence_panel.py:945-950``.
    """
    from backend.services.annotations import shift_annotations

    kept, dropped = shift_annotations(
        seq.annotations, start=start, end=end, inserted=len(inserted)
    )
    new_seq = seq.seq[:start] + inserted + seq.seq[end:]
    updated = seq.model_copy(
        update={"seq": new_seq, "annotations": kept, "length": len(new_seq)}
    )
    return updated, dropped


def flipped_sequence(seq):
    """Reverse-complement a record the way ``_commit_edit(flip=True)`` does."""
    from backend.services.annotations import flip_annotations

    rc = str(Seq(seq.seq).reverse_complement())
    return seq.model_copy(
        update={
            "seq": rc,
            "annotations": flip_annotations(seq.annotations, len(seq.seq)),
            "length": len(rc),
        }
    )


# ── a small GenBank record over real bases, for locations no real record has ──

def genbank_with(locations: list[tuple[str, str, str]], *, length: int = 300) -> str:
    """A minimal GenBank record carrying `locations` over real AF071878 bases.

    `locations` is a list of (feature_key, location_string, label). The bases
    are the first `length` of the pinned AF071878.1 record, so the input is real
    sequence; only the feature table is authored, because no public record in
    the corpus carries an ``order()`` or a remote reference.
    """
    seq = str(bio_record(BFDV_GB).seq)[:length]
    rows = []
    for i in range(0, length, 60):
        chunk = seq[i:i + 60]
        blocks = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        rows.append(f"{i + 1:>9} {blocks}")
    feats = "".join(
        f"     {key:<16}{loc}\n                     /label=\"{label}\"\n"
        for key, loc, label in locations
    )
    return (
        f"LOCUS       COORDPROBE{length:>18} bp    DNA     linear   SYN 01-JAN-2024\n"
        "DEFINITION  gx coords probe over real AF071878.1 bases.\n"
        "ACCESSION   COORDPROBE\n"
        "FEATURES             Location/Qualifiers\n"
        f"     source          1..{length}\n"
        f"{feats}"
        "ORIGIN\n" + "\n".join(rows) + "\n//\n"
    )


# ── fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def jcv_text() -> str:
    """The pinned NC_001699.1 GenBank bytes (see gx/corpus/MANIFEST.json)."""
    return JCV.read_text()


@pytest.fixture(scope="session")
def jcv_bio(jcv_text) -> SeqRecord:
    rec = SeqIO.read(StringIO(jcv_text), "genbank")
    assert len(rec.seq) == 5130 and len(rec.features) == 28
    return rec


@pytest.fixture
def jcv(jcv_text):
    """The subject's own parse of NC_001699.1 — a fresh model per test."""
    from backend.formats.genbank import read_genbank

    seqs = read_genbank(jcv_text)
    assert len(seqs) == 1
    return seqs[0]


@pytest.fixture(scope="session")
def bfdv_gb_text() -> str:
    return BFDV_GB.read_text()


@pytest.fixture(scope="session")
def bfdv_embl_text() -> str:
    return BFDV_EMBL.read_text()
