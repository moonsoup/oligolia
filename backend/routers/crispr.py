"""CRISPR guide RNA design endpoint."""

import re
from fastapi import APIRouter, HTTPException
from ..models.crispr import CRISPRDesignRequest, CRISPRDesignResponse, GuideRNA, CasType
from ..crispr_offtarget import scan_off_targets

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


#: The guide length each nuclease is normally used with, applied when the caller
#: does not ask for a specific one.
CANONICAL_GUIDE_LENGTH = {
    CasType.CAS9: 20,
    CasType.CAS9_HF: 20,
    CasType.CAS12A: 23,
    CasType.CAS13: 22,
}

#: How many bases of PAM have to fit alongside the guide. Cas13 targets RNA and
#: has no PAM requirement.
PAM_LENGTH = {
    CasType.CAS9: 3,       # NGG, 3' of the protospacer
    CasType.CAS9_HF: 3,
    CasType.CAS12A: 4,     # TTTV, 5' of the protospacer
    CasType.CAS13: 0,
}


@router.post("/design", response_model=CRISPRDesignResponse)
def design_guides(req: CRISPRDesignRequest) -> CRISPRDesignResponse:
    target = req.target_sequence.upper().replace(" ", "").replace("\n", "")
    if not all(c in "ACGTN" for c in target):
        raise HTTPException(400, "Target sequence must be DNA (ACGTN only)")
    cas = req.cas_type

    # `guide_length` was ignored by every branch — Cas9 hardcoded 20, Cas12a 23,
    # Cas13 22 (#64.3). It is honoured now when the caller actually sets it, and
    # otherwise each nuclease keeps its canonical length, so existing callers that
    # never passed it (and relied on 23 for Cas12a) are unaffected by the default
    # of 20 on the request model.
    guide_length = (
        req.guide_length
        if "guide_length" in req.model_fields_set
        else CANONICAL_GUIDE_LENGTH[cas]
    )

    # The PAM is part of what has to fit, and it differs per nuclease: NGG sits 3'
    # of the protospacer, TTTV sits 5' of it, and Cas13 has none. The old guard
    # always demanded guide_length + 3, which rejected a legitimate 22 nt Cas13
    # target as "too short".
    needed = guide_length + PAM_LENGTH[cas]
    if len(target) < needed:
        raise HTTPException(400, f"Target too short (need ≥{needed} nt for {cas.value})")

    guides: list[GuideRNA] = []

    if cas in (CasType.CAS9, CasType.CAS9_HF):
        # Forward strand: guide+PAM = 20 nt protospacer + N + GG.
        # The N BELONGS TO THE PAM. `(?=(.{20})GG)` captured it as the guide's
        # last base, so every guide was shifted one base along the target and the
        # real protospacer was never returned -- with the mismatch landing in the
        # PAM-proximal seed, where SpCas9 is least tolerant (#53).
        for m in re.finditer(rf"(?=(.{{{guide_length}}}).GG)", target):
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
        for m in re.finditer(rf"(?=(.{{{guide_length}}}).GG)", rc_target):
            guide_seq = m.group(1)
            # rc index s spans rc[s : s+L], which is target[len-s-L : len-s].
            pos = len(target) - m.start() - guide_length
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
        # PAM is 5'-TTTV-3' immediately followed by the guide. Both strands are
        # scanned: only the + strand was, so a target whose TTTV sites all sit on
        # the minus strand returned zero candidates rather than the real ones
        # (#64.2).
        pam_re = rf"(?=TTT[ACG](.{{{guide_length}}}))"
        for m in re.finditer(pam_re, target):
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
        rc_target = _reverse_complement(target)
        for m in re.finditer(pam_re, rc_target):
            guide_seq = m.group(1)
            # The guide occupies rc[s+4 : s+4+L]; map that back to the + strand.
            pos = len(target) - (m.start() + 4) - guide_length
            gc = _gc_content(guide_seq)
            guides.append(GuideRNA(
                sequence=guide_seq,
                pam="TTTV",
                position=pos,
                strand="-",
                gc_content=round(gc, 1),
                on_target_score=round(0.5 + (0.1 if 40 <= gc <= 70 else -0.1), 3),
            ))

    elif cas == CasType.CAS13:
        # Cas13 targets RNA, so every window of the given sense strand is a
        # candidate and there is no second strand to scan — the crRNA is
        # complementary to this transcript, not to its genomic reverse complement.
        # `range(0, len - L)` dropped the last window (#64.1).
        for i in range(0, len(target) - guide_length + 1):
            guide_seq = target[i:i + guide_length]
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

    if req.check_off_targets and cas != CasType.CAS13:
        # Off-target scanning applies to PAM-directed DNA nucleases. Cas13
        # targets RNA and has no genomic off-target model here, so it is skipped.
        cas_family = "cas12a" if cas == CasType.CAS12A else "cas9"
        references = [target, *req.reference_sequences]
        for g in guides:
            result = scan_off_targets(
                g.sequence,
                references,
                cas_family=cas_family,
                max_mismatches=req.max_mismatches,
            )
            g.off_target_count = result.total
            g.off_target_summary = result.summary
            g.specificity_score = result.specificity_score

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
