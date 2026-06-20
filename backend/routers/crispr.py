"""CRISPR guide RNA design endpoint."""

import re
from fastapi import APIRouter, HTTPException
from ..models.crispr import CRISPRDesignRequest, CRISPRDesignResponse, GuideRNA, CasType

router = APIRouter(prefix="/crispr", tags=["crispr"])


def _gc_content(seq: str) -> float:
    upper = seq.upper()
    return (upper.count("G") + upper.count("C")) / len(seq) * 100 if seq else 0.0


def _doench_rule_set1_score(guide: str) -> float:
    """
    Approximation of Doench et al. 2014 Rule Set 1 on-target scoring.
    Uses position-specific nucleotide preferences for 20-nt guides.
    Real RS1 uses logistic regression; this is a simplified heuristic.
    """
    guide = guide.upper()
    if len(guide) < 20:
        return 0.5
    score = 0.5
    # Preferred nucleotides at key positions (1-indexed as in the original paper)
    prefs = {
        3: {"A": 0.03, "T": 0.02},
        4: {"C": 0.03, "A": 0.02},
        10: {"G": 0.03},
        12: {"A": 0.03},
        13: {"G": 0.03},
        20: {"G": 0.02, "A": 0.02},
    }
    for pos, nucleotide_scores in prefs.items():
        nt = guide[pos - 1]
        score += nucleotide_scores.get(nt, 0)
    # GC content penalty (prefer 40-60%)
    gc = _gc_content(guide)
    if 40 <= gc <= 60:
        score += 0.1
    elif gc < 20 or gc > 80:
        score -= 0.2
    # Avoid poly-T stretches (reduces transcription)
    if "TTTT" in guide:
        score -= 0.15
    return round(min(max(score, 0.0), 1.0), 3)


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


@router.post("/design", response_model=CRISPRDesignResponse)
def design_guides(req: CRISPRDesignRequest) -> CRISPRDesignResponse:
    target = req.target_sequence.upper().replace(" ", "").replace("\n", "")
    if not all(c in "ACGTN" for c in target):
        raise HTTPException(400, "Target sequence must be DNA (ACGTN only)")
    if len(target) < req.guide_length + 3:
        raise HTTPException(400, f"Target too short (need ≥{req.guide_length + 3} nt)")

    cas = req.cas_type
    guides: list[GuideRNA] = []

    if cas in (CasType.CAS9, CasType.CAS9_HF):
        # Forward strand: guide+PAM = 20nt + NGG
        for m in re.finditer(r"(?=(.{20})GG)", target):
            guide_seq = m.group(1)
            pos = m.start()
            gc = _gc_content(guide_seq)
            score = _doench_rule_set1_score(guide_seq)
            guides.append(GuideRNA(
                sequence=guide_seq,
                pam="NGG",
                position=pos,
                strand="+",
                gc_content=round(gc, 1),
                on_target_score=score,
            ))
        # Reverse strand
        rc_target = _reverse_complement(target)
        for m in re.finditer(r"(?=(.{20})GG)", rc_target):
            guide_seq = m.group(1)
            pos = len(target) - m.start() - 20
            gc = _gc_content(guide_seq)
            score = _doench_rule_set1_score(guide_seq)
            guides.append(GuideRNA(
                sequence=guide_seq,
                pam="NGG",
                position=pos,
                strand="-",
                gc_content=round(gc, 1),
                on_target_score=score,
            ))

    elif cas == CasType.CAS12A:
        # PAM is 5'-TTTV-3' followed by 23 nt guide
        for m in re.finditer(r"(?=TTT[ACG](.{23}))", target):
            guide_seq = m.group(1)
            pos = m.start() + 4
            gc = _gc_content(guide_seq)
            guides.append(GuideRNA(
                sequence=guide_seq,
                pam="TTTV",
                position=pos,
                strand="+",
                gc_content=round(gc, 1),
                on_target_score=round(0.5 + (0.1 if 40 <= gc <= 70 else -0.1), 3),
            ))

    elif cas == CasType.CAS13:
        # CAS13 targets RNA; every 22 nt window is a candidate
        for i in range(0, len(target) - 22, 1):
            guide_seq = target[i:i+22]
            gc = _gc_content(guide_seq)
            guides.append(GuideRNA(
                sequence=guide_seq,
                pam="none",
                position=i,
                strand="+",
                gc_content=round(gc, 1),
                on_target_score=round(0.5 + (0.1 if 40 <= gc <= 60 else 0.0), 3),
            ))

    total = len(guides)
    # Sort by on-target score descending, then GC proximity to 50%
    guides.sort(key=lambda g: (-(g.on_target_score or 0), abs((g.gc_content or 0) - 50)))
    guides = guides[:req.max_guides]

    return CRISPRDesignResponse(
        target_length=len(target),
        cas_type=req.cas_type,
        guides=guides,
        total_candidates=total,
    )


@router.post("/score_guide")
def score_guide(guide_sequence: str) -> dict:
    """Score a single 20-nt guide RNA for on-target efficiency."""
    if len(guide_sequence) < 17 or len(guide_sequence) > 24:
        raise HTTPException(400, "Guide must be 17-24 nt")
    score = _doench_rule_set1_score(guide_sequence)
    gc = _gc_content(guide_sequence)
    issues = []
    if gc < 30:
        issues.append("Low GC content (<30%)")
    if gc > 75:
        issues.append("High GC content (>75%)")
    if "TTTT" in guide_sequence.upper():
        issues.append("Poly-T run (reduces RNA Pol III transcription)")
    if guide_sequence.upper().startswith("T"):
        issues.append("5' T may reduce efficiency with U6 promoter")
    return {
        "guide": guide_sequence,
        "on_target_score": score,
        "gc_content": round(gc, 1),
        "issues": issues,
        "recommendation": "Good" if score >= 0.6 else "Fair" if score >= 0.4 else "Poor",
    }
