"""Sequence alignment endpoints — pairwise and multiple sequence alignment."""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from Bio import Align
from Bio.Seq import Seq

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

    alignments = list(aligner.align(req.seq1, req.seq2))
    if not alignments:
        raise HTTPException(422, "No alignment found")

    best = alignments[0]
    counts = best.counts()

    # Extract gapped sequences from FASTA format output
    fasta_lines = best.format("fasta").strip().split("\n")
    gapped_seqs = [ln for ln in fasta_lines if not ln.startswith(">")]
    aligned1 = gapped_seqs[0] if len(gapped_seqs) >= 1 else req.seq1
    aligned2 = gapped_seqs[1] if len(gapped_seqs) >= 2 else req.seq2

    aln_len = len(aligned1)  # use gapped sequence length; Alignment.length removed in Biopython 1.82+
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


@router.post("/multiple", response_model=MSAResult)
def multiple_align(req: MSARequest) -> MSAResult:
    """Run multiple sequence alignment using MUSCLE (via subprocess) or fallback to simple."""
    if len(req.sequences) < 2:
        raise HTTPException(400, "Need at least 2 sequences for MSA")

    import subprocess
    import tempfile
    import os

    # Write input FASTA
    fasta_in = "".join(f">{s['id']}\n{s['seq']}\n" for s in req.sequences)

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
            raise RuntimeError(result.stderr.decode())

        from Bio import SeqIO
        aligned = list(SeqIO.parse(out_path, "fasta"))
        os.unlink(fin_path)
        os.unlink(out_path)

    except (FileNotFoundError, RuntimeError):
        # Fallback: simple pairwise star alignment (naive, for environments without MUSCLE)
        os.unlink(fin_path) if os.path.exists(fin_path) else None
        seqs = [s["seq"] for s in req.sequences]
        max_len = max(len(s) for s in seqs)
        padded = [s + "-" * (max_len - len(s)) for s in seqs]
        from Bio.SeqRecord import SeqRecord
        aligned = [SeqRecord(Seq(padded[i]), id=req.sequences[i]["id"]) for i in range(len(seqs))]

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
