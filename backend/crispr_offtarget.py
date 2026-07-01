"""Off-target prediction for CRISPR guide RNAs.

Scans one or more reference sequences for near-matches to a guide protospacer
adjacent to a valid PAM, buckets hits by mismatch count, and computes an
MIT-style specificity score (Hsu et al. 2013). This is an offline, sequence-based
approximation — it scans whatever reference context the caller provides (the
target itself by default, plus any additional loci), not a whole genome.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# Hsu et al. 2013 position-specific mismatch weights for SpCas9.
# Index 0 == protospacer position 1 (PAM-distal 5' end);
# index 19 == position 20 (PAM-proximal, the intolerant seed).
MIT_WEIGHTS: list[float] = [
    0.000, 0.000, 0.014, 0.000, 0.000,
    0.395, 0.317, 0.000, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732,
    0.828, 0.615, 0.804, 0.685, 0.583,
]

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class OffTargetResult:
    """Off-target summary for a single guide."""

    summary: dict[str, int] = field(default_factory=dict)  # {"0": n, "1": n, ...}
    total: int = 0
    specificity_score: float = 100.0


def _mismatch_positions(guide: str, candidate: str) -> list[int]:
    """1-indexed positions (5'→3') where guide and candidate differ."""
    return [i + 1 for i, (a, b) in enumerate(zip(guide, candidate)) if a != b]


def _mit_hit_score(mm_positions: list[int]) -> float:
    """MIT single-site off-target score in [0, 1] (1.0 == identical site)."""
    if not mm_positions:
        return 1.0
    n = len(mm_positions)
    # Term 1: product of per-position tolerance for each mismatch.
    t1 = 1.0
    for pos in mm_positions:
        idx = min(pos - 1, len(MIT_WEIGHTS) - 1)
        t1 *= 1.0 - MIT_WEIGHTS[idx]
    # Term 2: mean pairwise distance between mismatches (clustered = worse).
    if n > 1:
        pairs = [
            abs(a - b)
            for i, a in enumerate(mm_positions)
            for b in mm_positions[i + 1:]
        ]
        d = sum(pairs) / len(pairs)
    else:
        d = 19.0
    t2 = 1.0 / (((19.0 - d) / 19.0) * 4.0 + 1.0)
    # Term 3: penalty for number of mismatches.
    t3 = 1.0 / (n * n)
    return t1 * t2 * t3


def _pam_positions_cas9(reference: str, guide_len: int) -> list[int]:
    """Start indices of guide-length windows immediately followed by a GG PAM.

    Matches the convention used by the guide designer (``routers/crispr.py``):
    the protospacer is a ``guide_len``-mer whose next two bases are ``GG``.
    """
    positions = []
    for i in range(len(reference) - guide_len - 2 + 1):
        if reference[i + guide_len] == "G" and reference[i + guide_len + 1] == "G":
            positions.append(i)
    return positions


def _pam_positions_cas12a(reference: str, guide_len: int) -> list[int]:
    """Start indices of the guide window for a 5' TTTV PAM."""
    positions = []
    for i in range(len(reference) - guide_len - 4 + 1):
        if (
            reference[i] == "T"
            and reference[i + 1] == "T"
            and reference[i + 2] == "T"
            and reference[i + 3] in "ACG"
        ):
            positions.append(i + 4)  # guide starts right after TTTV
    return positions


def scan_off_targets(
    guide: str,
    references: list[str],
    *,
    cas_family: str = "cas9",
    max_mismatches: int = 3,
) -> OffTargetResult:
    """Scan references for off-target sites of ``guide``.

    ``references`` must include the guide's own origin sequence exactly once;
    the single guaranteed on-target (0-mismatch) hit is subtracted so the
    result reflects off-targets only. Both strands are searched.
    """
    guide = guide.upper()
    guide_len = len(guide)
    summary = {str(k): 0 for k in range(max_mismatches + 1)}

    for ref in references:
        ref = ref.upper().replace(" ", "").replace("\n", "")
        for strand_seq in (ref, _reverse_complement(ref)):
            if cas_family == "cas12a":
                starts = _pam_positions_cas12a(strand_seq, guide_len)
            else:
                starts = _pam_positions_cas9(strand_seq, guide_len)
            for start in starts:
                candidate = strand_seq[start:start + guide_len]
                if len(candidate) != guide_len:
                    continue
                mm = _mismatch_positions(guide, candidate)
                if len(mm) <= max_mismatches:
                    summary[str(len(mm))] += 1

    # Remove the guide's own on-target: it matches its origin with 0 mismatches.
    if summary.get("0", 0) > 0:
        summary["0"] -= 1

    # Aggregate MIT specificity from the surviving off-target hits.
    hit_score_sum = 0.0
    for n_mm_str, count in summary.items():
        if count <= 0:
            continue
        # Approximate: treat all hits in a bucket as the worst case for that
        # bucket (mismatches in the most-tolerant PAM-distal positions).
        positions = list(range(1, int(n_mm_str) + 1)) or []
        per_hit = _mit_hit_score(positions)
        hit_score_sum += per_hit * count

    specificity = 100.0 / (1.0 + hit_score_sum)

    total = sum(summary.values())
    return OffTargetResult(
        summary=summary,
        total=total,
        specificity_score=round(specificity, 1),
    )
