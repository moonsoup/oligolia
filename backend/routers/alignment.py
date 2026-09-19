"""Sequence alignment endpoints — pairwise and multiple sequence alignment."""

import shutil

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


#: The aligner used unless a caller asks for another, and the one the GUI drives.
#: Named here so the panel can describe it without hard-coding its name (#83).
DEFAULT_ALGORITHM = "muscle"


class MSARequest(BaseModel):
    sequences: list[dict]  # [{"id": str, "seq": str}]
    algorithm: str = DEFAULT_ALGORITHM  # muscle | clustalw


class MSAResult(BaseModel):
    aligned: list[dict]  # [{"id": str, "aligned_seq": str}]
    consensus: str
    identity_matrix: list[list[float]]


class AlignerInfo(BaseModel):
    """Whether one external aligner can be run, and what to do if it cannot."""

    name: str
    available: bool
    path: str | None  # resolved executable, or None when nothing was found
    hint: str  # the same install hint the 503 quotes
    url: str  # where to get it


class AlignersResult(BaseModel):
    aligners: list[AlignerInfo]
    any_available: bool


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


#: Homepage per aligner, quoted inside the hint below and offered on its own so a
#: GUI can turn it into a link without parsing the sentence back apart.
_ALIGNER_URLS = {
    "muscle": "https://drive5.com/muscle",
    "clustalw": "http://www.clustal.org/clustal2",
}

#: Where to get each aligner, so a 503 can tell the user what to do.
_ALIGNER_HELP = {
    "muscle": f"MUSCLE v5 (`brew install muscle`, or {_ALIGNER_URLS['muscle']})",
    "clustalw": f"ClustalW (`brew install clustal-w`, or {_ALIGNER_URLS['clustalw']})",
}


def aligner_status(name: str) -> AlignerInfo:
    """Is `name` runnable right now, and what to install if it is not.

    The one place that answers the question. `/multiple` asks it before shelling
    out and `/aligners` (and through it the MSA panel) asks it up front, so the
    tab cannot say "ready" about an aligner the run path would 503 on, and cannot
    advise an install the error names differently (#83).
    """
    path = shutil.which(name)
    return AlignerInfo(
        name=name,
        available=path is not None,
        path=path,
        hint=_ALIGNER_HELP[name],
        url=_ALIGNER_URLS[name],
    )


def aligner_statuses() -> list[AlignerInfo]:
    """`aligner_status` for every aligner `/multiple` knows how to drive."""
    return [aligner_status(name) for name in _ALIGNER_HELP]


def _not_installed(algorithm: str) -> HTTPException:
    """The 503 for a missing aligner — worded once, raised from both paths."""
    return HTTPException(
        503,
        f"{algorithm} is not installed, so these sequences cannot be aligned. "
        f"Install {_ALIGNER_HELP[algorithm]} and try again. "
        "Nothing is returned rather than an approximation, because a padded "
        "copy of the input is indistinguishable from a real alignment (#58).",
    )


@router.get("/aligners", response_model=AlignersResult)
def list_aligners() -> AlignersResult:
    """Report aligner availability without running an alignment.

    Additive: `/multiple` behaves exactly as it did, including its 503. This
    exists so the MSA tab can disable itself and name the install *before* the
    user pastes sequences and presses Run, rather than after (#83, from #76's
    option C).
    """
    statuses = aligner_statuses()
    return AlignersResult(
        aligners=statuses,
        any_available=any(a.available for a in statuses),
    )


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
    if not aligner_status(req.algorithm).available:
        # Same answer the panel got up front, same 503 the FileNotFoundError arm
        # below raises — asking first only means it is raised without writing a
        # temp file first. The arm below stays as the backstop for an executable
        # that exists on PATH but cannot actually be exec'd.
        raise _not_installed(req.algorithm)

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
        raise _not_installed(req.algorithm) from None
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
