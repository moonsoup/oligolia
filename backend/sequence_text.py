"""One definition of what a pasted sequence *is* (#86).

Every endpoint and GUI call site that analyses the same template box has to
agree on which characters are the molecule and which are layout, or the same
paste gets two different sets of coordinates — the Restriction tab's positions
and the Digest tab's cut/start/end columns come from different endpoints.
"""

from __future__ import annotations

import string

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
    normalises, and rejecting it is a separate concern (#87).
    """
    return text.translate(_STRIP_LAYOUT).upper()
