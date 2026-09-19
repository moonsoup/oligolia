"""One definition of what a pasted sequence *is* (#86).

Every endpoint and GUI call site that analyses the same template box has to
agree on which characters are the molecule and which are layout, or the same
paste gets two different sets of coordinates — the Restriction tab's positions
and the Digest tab's cut/start/end columns come from different endpoints.
"""

from __future__ import annotations

import string

from fastapi import HTTPException

# Exactly the characters Bio.Restriction's FormattedSeq removes before it
# searches — Biopython 1.85, Bio/Restriction/Restriction.py:
#     _remove_chars = string.whitespace.encode() + string.digits.encode()
# /primers/digest's cut positions come back from that search, so this is the
# coordinate space those numbers already live in. Stripping anything less (the
# old ``.replace(" ", "").replace("\n", "")`` left \r, \t and the ORIGIN line
# numbers behind) means slicing a string Biopython never saw.
_LAYOUT_CHARS = string.whitespace + string.digits
_STRIP_LAYOUT = str.maketrans("", "", _LAYOUT_CHARS)


def normalize_template(text: str) -> str:
    """Uppercase ``text`` with sequence-layout characters removed.

    A GenBank ORIGIN block (numbered, spaced, 60-column), a CRLF-wrapped FASTA
    body and the bare sequence all denote the same molecule, and this maps all
    three onto it. Note that a line break falling inside a recognition site is
    the same case: the site is only found once the break is gone.

    Characters that are *not* layout — a leftover FASTA header, an ``X`` — are
    deliberately left in place. That is an input Biopython refuses rather than
    normalises, so it is not this function's job to hide it; use
    :func:`validate_template` to turn it into a refusal (#87).
    """
    return text.translate(_STRIP_LAYOUT).upper()


# Exactly the letters Bio.Restriction's FormattedSeq accepts once layout is
# gone — Biopython 1.85, Bio/Restriction/Restriction.py:
#     for c in b"ABCDGHKMNRSTVWY":  # Only allow IUPAC letters
# Anything else and it raises TypeError("Invalid character found in ...") rather
# than guessing, so a template containing one is not an analysable sequence.
#
# This is Bio.Data.IUPACData.ambiguous_dna_values' alphabet minus ``X``:
# IUPACData lists "X": "GATC" as a legacy any-base code, but FormattedSeq's
# table does not, and the authority here has to be the code that does the
# searching — validating against ``X`` would pass an input that then crashes.
IUPAC_NUCLEOTIDES = frozenset("ABCDGHKMNRSTVWY")


def first_non_iupac(seq: str) -> tuple[int, str] | None:
    """``(0-based index, character)`` of the first non-IUPAC letter, or ``None``.

    ``seq`` is expected to be already normalised by :func:`normalize_template`,
    so layout characters have been removed and lowercase upper-cased; the index
    is therefore into the molecule as the app reads it.
    """
    for i, char in enumerate(seq):
        if char not in IUPAC_NUCLEOTIDES:
            return i, char
    return None


def non_iupac_detail(index: int, char: str, label: str = "Template") -> str:
    """The 400 message for a template that is not a nucleotide sequence."""
    detail = (
        f"{label} is not a nucleotide sequence: invalid character {char!r} at "
        f"position {index + 1}. Expected IUPAC nucleotide codes "
        f"(A, C, G, T plus ambiguity codes {', '.join(sorted(IUPAC_NUCLEOTIDES - set('ACGT')))})."
    )
    if char == ">":
        detail += " If you pasted a FASTA record, delete its '>' description line."
    return detail


def validate_template(text: str, label: str = "Template") -> str:
    """Normalise ``text`` and refuse it if it is not a nucleotide sequence.

    One place where both ``/primers/restriction_sites`` and ``/primers/digest``
    decide whether their input is a sequence at all, so a paste that still has
    its FASTA header gets the same HTTP 400 from each instead of one endpoint
    raising an unhandled ``TypeError`` out of Biopython (500) and the other
    returning 200 with positions counted over the header text (#87).

    Raises ``fastapi.HTTPException`` (400), matching how a bad enzyme name is
    already rejected, so the GUI's existing error paths surface a message about
    the user's sequence rather than a Biopython internal.
    """
    seq = normalize_template(text)
    bad = first_non_iupac(seq)
    if bad is not None:
        raise HTTPException(400, non_iupac_detail(*bad, label=label))
    return seq
