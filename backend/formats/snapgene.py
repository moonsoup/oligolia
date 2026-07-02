"""SnapGene ``.dna`` binary format import (issue #41).

SnapGene's native format is binary. Biopython's ``SnapGene`` SeqIO parser
already handles it (no new dependency), and its SeqRecords carry the same
shape as GenBank — sequence, features with locations, and a ``topology``
annotation — so conversion reuses the GenBank ``_to_sequence`` logic and
therefore inherits circular/linear topology handling (issue #40).

Import only; SnapGene's own GenBank export is the practical path back out.
"""

from __future__ import annotations

from io import BytesIO
from typing import BinaryIO

from Bio import SeqIO

from ..models.sequence import Sequence
from .genbank import _to_sequence


def read_snapgene(source: str | bytes | BinaryIO) -> list[Sequence]:
    """Parse a SnapGene ``.dna`` file into Sequence objects.

    ``source`` may be a filesystem path, a binary file handle, or raw bytes.
    A SnapGene file holds a single sequence, but a list is returned for
    consistency with the other format readers.
    """
    if isinstance(source, (bytes, bytearray)):
        handle: object = BytesIO(source)
    else:
        handle = source  # path string or binary file handle
    records = list(SeqIO.parse(handle, "snapgene"))
    return [_to_sequence(r) for r in records]
