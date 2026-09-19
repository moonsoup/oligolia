"""Independent oracles for the gx ORF/translation round.

Imported as ``from orf_oracle import ...`` (never ``from conftest import ...``:
that is what forced earlier rounds into separate processes). ``conftest.py`` in
this directory holds pytest fixtures only.

Nothing here reads the subject. The three sources of truth are:

* ``gx/corpus/ncbi/gc.prt`` -- NCBI's own genetic-code tables, the human-readable
  authority for what table 1/2/11 map and which codons each calls a start.
* ``Bio.Data.CodonTable`` / ``Bio.Seq.translate`` (Biopython 1.85, pinned at
  ``gx/corpus/biopython-1.85/Bio/Data/CodonTable.py`` and ``.../Bio/Seq.py``) --
  the machine-readable form, cross-checked against gc.prt by
  :func:`gc_prt_table` at import of the tests that use it.
* the ``/translation`` qualifiers NCBI ships inside the pinned records, which are
  recorded literals: nothing in this repo produced them.
"""

from __future__ import annotations

import re
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq

CORPUS = Path(__file__).resolve().parent.parent / "corpus"
GC_PRT = CORPUS / "ncbi" / "gc.prt"

#: Codon order of the gc.prt ``ncbieaa``/``sncbieaa`` strings, spelled out by the
#: three ``-- Base1/2/3`` comment lines that follow every table in that file.
_GC_PRT_BASES = "TCAG"
GC_PRT_CODON_ORDER: list[str] = [
    b1 + b2 + b3 for b1 in _GC_PRT_BASES for b2 in _GC_PRT_BASES for b3 in _GC_PRT_BASES
]


def gc_prt_table(table_id: int) -> tuple[str, str]:
    """``(ncbieaa, sncbieaa)`` for ``table_id``, read out of the pinned gc.prt.

    ``ncbieaa`` is the amino acid of each of the 64 codons in
    :data:`GC_PRT_CODON_ORDER`, ``*`` for a stop; ``sncbieaa`` marks a start
    codon with ``M`` and everything else with ``-`` (stops stay ``*``).
    """
    text = GC_PRT.read_text()
    m = re.search(
        r"id %d ,\s*ncbieaa\s+\"([^\"]*)\"\s*,\s*sncbieaa\s+\"([^\"]*)\"" % table_id,
        text,
    )
    if m is None:  # pragma: no cover - a corpus regression, not a subject one
        raise AssertionError(f"no table id {table_id} in {GC_PRT}")
    return m.group(1), m.group(2)


def ncbi_amino_acids(table_id: int) -> dict[str, str]:
    """``{codon: aa}`` over all 64 codons straight from gc.prt, ``*`` for stop."""
    ncbieaa, _ = gc_prt_table(table_id)
    return dict(zip(GC_PRT_CODON_ORDER, ncbieaa, strict=True))


def ncbi_start_codons(table_id: int) -> set[str]:
    """The start codons gc.prt declares for ``table_id``.

    Table 1 is ``{TTG, CTG, ATG}``; table 11 adds ``GTG``, ``ATT``, ``ATC``,
    ``ATA``. This is the set an ``include_all_starts`` option has to be drawn
    from -- a start codon is a property of a table, not a free-floating list.
    """
    _, sncbieaa = gc_prt_table(table_id)
    return {c for c, s in zip(GC_PRT_CODON_ORDER, sncbieaa, strict=True) if s == "M"}


def biopython_matches_gc_prt(table_id: int) -> bool:
    """True when Biopython's table ``table_id`` renders back to gc.prt exactly.

    The two authorities are independent (one is NCBI's ASN.1 text, the other is
    Biopython's hand-maintained Python), so agreement between them is what makes
    either usable as an oracle for the subject.
    """
    t = CodonTable.unambiguous_dna_by_id[table_id]
    ncbieaa = "".join(
        "*" if c in t.stop_codons else t.forward_table[c] for c in GC_PRT_CODON_ORDER
    )
    sncbieaa = "".join(
        "*" if c in t.stop_codons else ("M" if c in t.start_codons else "-")
        for c in GC_PRT_CODON_ORDER
    )
    return (ncbieaa, sncbieaa) == gc_prt_table(table_id)


# ── Re-extraction by index ────────────────────────────────────────────────────

def orf_bases(sequence: str, start: int, end: int, frame: int) -> str:
    """The bases an ORF reported as ``[start, end)`` on ``frame`` covers.

    ``start``/``end`` are 0-based half-open on the TOP strand for both strands --
    that is what ``backend/routers/analysis.ORF.start`` documents ("0-based
    position on input strand"). A negative frame is read on the reverse
    complement of that same top-strand slice, which is the only mapping under
    which the reported ``start_codon`` can be at the 5' end of the ORF.
    """
    piece = sequence[start:end]
    return str(Seq(piece).reverse_complement()) if frame < 0 else piece


def translate_codons(bases: str, table_id: int = 1) -> str:
    """Codon-by-codon translation of ``bases`` under ``table_id``, stops as ``*``.

    Any trailing 1-2 bases are NOT translated and NOT silently absorbed: the
    caller is handed ``len(bases) % 3`` separately by :func:`partial_tail`.
    Deliberately does not call ``Seq.translate`` on a ragged string, so the
    BiopythonWarning this round is about can never be raised from the oracle.
    """
    aa = ncbi_amino_acids(table_id)
    whole = bases[: len(bases) // 3 * 3]
    return "".join(aa.get(whole[i:i + 3], "X") for i in range(0, len(whole), 3))


def partial_tail(bases: str) -> str:
    """The 1-2 trailing bases that are not a whole codon (``''`` when aligned)."""
    return bases[len(bases) // 3 * 3:]


def biopython_translation(bases: str, table_id: int = 1) -> str:
    """``Bio.Seq.translate`` over whole codons only, stops as ``*``.

    Trimmed before the call for the same reason as above -- the oracle must not
    depend on Biopython's own partial-codon behaviour while the round is asking
    whether the subject handles it.
    """
    whole = bases[: len(bases) // 3 * 3]
    return str(Seq(whole).translate(table=table_id))


def frame_substring(sequence: str, frame: int) -> str:
    """The bases reading frame ``frame`` (+/-1, +/-2, +/-3) translates.

    Positive frames read ``sequence[offset:]``; negative frames read
    ``reverse_complement(sequence)[offset:]``, offset ``abs(frame) - 1``. This is
    the frame decomposition ``find_orfs`` itself walks
    (``backend/routers/analysis.py:509-511``), stated independently of it.
    """
    offset = abs(frame) - 1
    strand = sequence if frame > 0 else str(Seq(sequence).reverse_complement())
    return strand[offset:]


def top_strand_span(seq_len: int, strand_start: int, strand_end: int, frame: int) -> tuple[int, int]:
    """Map a half-open span in frame-strand space back to the top strand.

    For a positive frame the span is already top-strand. For a negative frame a
    span ``[s, e)`` on the reverse complement of an ``n``-base molecule is
    ``[n - e, n - s)`` on the top strand -- the mapping is its own inverse, so
    the round trip is checkable.
    """
    if frame > 0:
        return strand_start, strand_end
    return seq_len - strand_end, seq_len - strand_start


# ── An ORF finder that is not the subject's ───────────────────────────────────

def find_orfs_oracle(
    sequence: str,
    *,
    min_length_aa: int = 30,
    table_id: int = 1,
    start_codons: set[str] | None = None,
) -> list[dict]:
    """Re-derive the six-frame ORFs of a LINEAR ``sequence`` from gc.prt alone.

    Written from the definition rather than from the subject: for each of the six
    frames, walk codons; on a start codon open an ORF; on a stop codon close it
    and record it if it is long enough; at the end of the frame close whatever is
    still open and mark it ``partial=True``.

    ``start`` / ``end`` are 0-based half-open TOP-STRAND coordinates and always
    delimit whole translated codons -- an ORF left open at the end of a frame
    stops at its last complete codon, so the 1-2 untranslated trailing bases fall
    OUTSIDE the span rather than inside it.

    ``protein`` is the literal translation under ``table_id``: the initiator is
    whatever that table says the codon is, with no initiator-methionine rewrite.
    A caller that wants the GenBank ``/translation`` convention applies
    :func:`as_initiator_met` to the result and says so.

    Circularity is deliberately absent. The circular case is tested by rotation
    invariance instead (``test_orf_circular.py``), which needs no second ORF
    finder: the same molecule read from a different origin must yield the same
    ORFs, and that is checkable by index against the unrotated run.
    """
    aa_of = ncbi_amino_acids(table_id)
    starts = ncbi_start_codons(table_id) if start_codons is None else set(start_codons)
    n = len(sequence)
    rc = str(Seq(sequence).reverse_complement())

    out: list[dict] = []
    for sign in (+1, -1):
        strand = sequence if sign > 0 else rc
        for offset in range(3):
            i = offset
            open_at: int | None = None
            residues: list[str] = []
            while i + 3 <= len(strand):
                codon = strand[i:i + 3]
                aa = aa_of.get(codon, "X")
                if open_at is None:
                    if codon in starts:
                        open_at, residues = i, [aa]
                elif aa == "*":
                    if len(residues) >= min_length_aa:
                        out.append(_record(sign, offset, open_at, i + 3, residues, n, False, strand))
                    open_at, residues = None, []
                else:
                    residues.append(aa)
                i += 3
            if open_at is not None and len(residues) >= min_length_aa:
                out.append(_record(sign, offset, open_at, i, residues, n, True, strand))
    out.sort(key=lambda o: (-o["length_aa"], o["frame"], o["start"]))
    return out


def _record(sign, offset, s, e, residues, n, partial, strand) -> dict:
    start, end = top_strand_span(n, s, e, sign)
    return {
        "frame": sign * (offset + 1),
        "start": start,
        "end": end,
        "length_nt": e - s,
        "length_aa": len(residues),
        "protein": "".join(residues),
        "start_codon": strand[s:s + 3],
        "partial": partial,
    }


def as_initiator_met(protein: str) -> str:
    """``protein`` with its first residue rendered ``M``.

    The GenBank ``/translation`` convention: whatever the initiator codon codes
    for elsewhere in the chain, the first residue of a CDS is written ``M``. See
    ``gx/corpus/inputs/NC_012920.1.gb``'s MT-ND2, whose ``ATT`` start is ``I`` in
    a literal table-2 translation and ``M`` in the record NCBI ships.
    """
    return "M" + protein[1:] if protein else protein
