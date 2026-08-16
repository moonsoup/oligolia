"""ESMFold public REST API client — https://esmatlas.com/about (no API key required).

Free public inference endpoint for single-chain protein structure prediction.
Last-resort fallback: only used when a sequence has no UniProt accession match
in AlphaFold DB (alphafold_db.py) — i.e. novel, synthetic, or in-app-mutated
sequences that can't already have a precomputed structure.

Practical length ceiling ~400 residues (server-enforced); longer sequences are
rejected client-side with a clear message rather than silently truncated.
"""

import httpx

BASE = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_LENGTH = 400


class ESMFoldClient:
    def __init__(self, timeout: int = 120) -> None:
        self._timeout = timeout

    def predict(self, sequence: str) -> str:
        """Fold a raw protein sequence. Returns PDB text with per-residue pLDDT
        confidence encoded in the B-factor column.

        Raises ValueError (not a network error) for sequences over MAX_LENGTH —
        the caller should surface this as a specific, actionable message rather
        than a generic failure.
        """
        seq = sequence.strip().upper()
        if len(seq) > MAX_LENGTH:
            raise ValueError(
                f"Sequence is {len(seq)} residues — ESMFold's free API supports "
                f"up to {MAX_LENGTH}. Try a shorter fragment, or look it up by "
                f"UniProt accession instead (AlphaFold DB has no length limit)."
            )
        with httpx.Client(timeout=self._timeout) as client:
            r = client.post(BASE, content=seq)
            r.raise_for_status()
            return r.text
