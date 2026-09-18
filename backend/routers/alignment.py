"""Sequence alignment endpoints — pairwise and multiple sequence alignment."""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from Bio import Align

router = APIRouter(prefix="/alignment", tags=["alignment"])


class PairwiseRequest(BaseModel):
    seq1: str
    seq2: str
    mode: str = "global"  # global | local
    match_score: float = 2.0
    mismatch_score: float = -1.0
    open_gap_score: float = -2.0
    extend_gap_score: float = -0.5


class PairwiseResult(BaseModel):
    score: float
    aligned_seq1: str
    aligned_seq2: str
    identity: float
    similarity: float
    gaps: int
    alignment_length: int


class MSARequest(BaseModel):
    sequences: list[dict]  # [{"id": str, "seq": str}]
    algorithm: str = "muscle"  # muscle | clustalw


class MSAResult(BaseModel):
    aligned: list[dict]  # [{"id": str, "aligned_seq": str}]
    consensus: str
    identity_matrix: list[list[float]]


@router.post("/pairwise", response_model=PairwiseResult)
def pairwise_align(req: PairwiseRequest) -> PairwiseResult:
    aligner = Align.PairwiseAligner()
    aligner.mode = req.mode
    aligner.match_score = req.match_score
    aligner.mismatch_score = req.mismatch_score
    aligner.open_gap_score = req.open_gap_score
    aligner.extend_gap_score = req.extend_gap_score

    # Take the first (optimal) alignment lazily. Never materialise the set and never
    # call len() on it: two unrelated 500-nt sequences already have more co-optimal
    # alignments than fit in an int64, so `list(...)` raised OverflowError/MemoryError
    # on ordinary input while only alignments[0] was ever used (#56).
    try:
        best = next(iter(aligner.align(req.seq1, req.seq2)))
    except StopIteration:
        raise HTTPException(422, "No alignment found") from None
    counts = best.counts()

    # Extract gapped sequences from FASTA format output
    fasta_lines = best.format("fasta").strip().split("\n")
    gapped_seqs = [ln for ln in fasta_lines if not ln.startswith(">")]
    aligned1 = gapped_seqs[0] if len(gapped_seqs) >= 1 else req.seq1
    aligned2 = gapped_seqs[1] if len(gapped_seqs) >= 2 else req.seq2

    aln_len = len(aligned1)  # gapped sequence length (equals best.length for these modes)
    identity = counts.identities / aln_len if aln_len else 0
    similarity = (counts.identities + counts.mismatches) / aln_len if aln_len else 0

    return PairwiseResult(
        score=best.score,
        aligned_seq1=aligned1,
        aligned_seq2=aligned2,
        identity=round(identity * 100, 2),
        similarity=round(similarity * 100, 2),
        gaps=counts.gaps,
        alignment_length=aln_len,
    )


#: Where to get each aligner, so a 503 can tell the user what to do.
_ALIGNER_HELP = {
    "muscle": "MUSCLE v5 (`brew install muscle`, or https://drive5.com/muscle)",
    "clustalw": "ClustalW (`brew install clustal-w`, or http://www.clustal.org/clustal2)",
}


@router.post("/multiple", response_model=MSAResult)
def multiple_align(req: MSARequest) -> MSAResult:
    """Run multiple sequence alignment with an external aligner.

    Refuses rather than approximating. The previous fallback right-padded the input
    with "-" and returned it as an alignment, complete with a consensus and an
    identity matrix — so two sequences differing by a one-base offset came back at
    0% identity while the UI said "Aligned N sequences" (#58). Neither aligner is
    bundled, so that was the normal path for users rather than an edge case.

    A real built-in aligner is wanted and is tracked separately; an approximation
    that cannot be told apart from a true alignment is worse than an error, which is
    the whole lesson of #58.
    """
    if len(req.sequences) < 2:
        raise HTTPException(400, "Need at least 2 sequences for MSA")
    if req.algorithm not in _ALIGNER_HELP:
        raise HTTPException(
            400, f"Unknown algorithm {req.algorithm!r}; expected one of {sorted(_ALIGNER_HELP)}"
        )

    import subprocess
    import tempfile
    import os

    # Write input FASTA
    fasta_in = "".join(f">{s['id']}\n{s['seq']}\n" for s in req.sequences)

    fin_path: str | None = None
    out_path: str | None = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fin:
            fin.write(fasta_in)
            fin_path = fin.name
        out_path = fin_path + ".aln"

        if req.algorithm == "muscle":
            result = subprocess.run(
                ["muscle", "-align", fin_path, "-output", out_path],
                capture_output=True, timeout=60,
            )
        else:  # clustalw
            result = subprocess.run(
                ["clustalw", "-INFILE=" + fin_path, "-OUTFILE=" + out_path, "-OUTPUT=FASTA"],
                capture_output=True, timeout=60,
            )

        if result.returncode != 0:
            raise HTTPException(
                502,
                f"{req.algorithm} ran but failed (exit {result.returncode}): "
                f"{result.stderr.decode(errors='replace').strip()[:500]}",
            )

        from Bio import SeqIO
        aligned = list(SeqIO.parse(out_path, "fasta"))
        if len(aligned) != len(req.sequences):
            raise HTTPException(
                502,
                f"{req.algorithm} returned {len(aligned)} sequences for "
                f"{len(req.sequences)} inputs — the alignment is incomplete",
            )

    except FileNotFoundError:
        raise HTTPException(
            503,
            f"{req.algorithm} is not installed, so these sequences cannot be aligned. "
            f"Install {_ALIGNER_HELP[req.algorithm]} and try again. "
            "Nothing is returned rather than an approximation, because a padded "
            "copy of the input is indistinguishable from a real alignment (#58).",
        ) from None
    except subprocess.TimeoutExpired:
        # Was uncaught, so a slow aligner produced a 500 with no explanation.
        raise HTTPException(
            504,
            f"{req.algorithm} did not finish within 60s for "
            f"{len(req.sequences)} sequences (longest {max(len(s['seq']) for s in req.sequences)} bp)",
        ) from None
    finally:
        # The old code leaked out_path on every fallback and fin_path on some.
        for p in (fin_path, out_path):
            if p and os.path.exists(p):
                os.unlink(p)

    aligned_out = [{"id": r.id, "aligned_seq": str(r.seq)} for r in aligned]

    # Consensus
    if aligned:
        length = len(aligned[0].seq)
        consensus = ""
        for i in range(length):
            col = [str(r.seq[i]).upper() for r in aligned]
            most = max(set(col), key=col.count)
            consensus += most if col.count(most) > len(col) / 2 else "N"
    else:
        consensus = ""

    # Pairwise identity matrix
    n = len(aligned)
    matrix: list[list[float]] = [[0.0] * n for _ in range(n)]
    for i in range(n):
        matrix[i][i] = 100.0
        for j in range(i + 1, n):
            s1, s2 = str(aligned[i].seq), str(aligned[j].seq)
            same = sum(a == b for a, b in zip(s1, s2) if a != "-" and b != "-")
            total = sum(1 for a, b in zip(s1, s2) if a != "-" or b != "-")
            pct = round(same / total * 100, 2) if total else 0.0
            matrix[i][j] = matrix[j][i] = pct

    return MSAResult(aligned=aligned_out, consensus=consensus, identity_matrix=matrix)
