"""INDEPENDENT oracles for the gx CRISPR round. Imported as ``crispr_oracle``.

Deliberately NOT ``conftest``: every gx round that did ``from conftest import ...``
is why rounds have to run in separate processes. ``gx/crispr/conftest.py`` holds
pytest fixtures and nothing else.

Nothing here calls ``backend.routers.crispr`` or ``backend.crispr_offtarget``.
Every expected value comes from one of four places, never from the subject:

* a regex scan of the raw sequence written from the PAM grammar, not from the
  subject's own patterns;
* re-derivation by index out of the raw string (``seq[i:j]``);
* Biopython 1.85 -- ``Bio.Seq.reverse_complement`` and
  ``Bio.SeqUtils.gc_fraction``, pinned at
  ``gx/corpus/biopython-1.85/Bio/Seq.py`` and
  ``gx/corpus/biopython-1.85/Bio/SeqUtils/__init__.py``;
* the two published scoring definitions the subject names, pinned at
  ``gx/corpus/crispr/crisporEffScores-v3.1.py`` (Doench et al. 2014 Rule Set 1)
  and ``gx/corpus/crispr/crispor-v3.1.py`` (Hsu et al. 2013 / MIT).

The two scoring tables below are re-typed rather than imported -- the pinned
CRISPOR files are reference bytes under a contested licence (see
gx/corpus/MANIFEST.json), not a dependency. ``_assert_against_pinned_bytes()``
runs at import time and fails the whole suite if a re-typed number has drifted
from the pinned copy, so the transcription cannot silently rot.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

from Bio.Seq import reverse_complement
from Bio.SeqUtils import gc_fraction

GX = Path(__file__).resolve().parent.parent
CORPUS = GX / "corpus"
INPUTS = CORPUS / "inputs"
REGISTER = json.loads((Path(__file__).resolve().parent / "register.json").read_text())


def requirement(rid: str) -> dict:
    """The pinned requirement, so tests read their numbers from the register."""
    for r in REGISTER:
        if r["id"] == rid:
            return r
    raise AssertionError(f"no requirement {rid!r} in gx/crispr/register.json")


# ── the PAM grammar, written from the authority and not from the subject ─────
#
# SpCas9 / SpCas9-HF1: protospacer, then NGG on its 3' side.
#   gx/corpus/crispr/crispor-v3.1.py#('NGG','NGG - Streptococcus Pyogenes')
#   Cong et al. 2013 (PMID 23287718, PMC3795411): "a 30-base pair (bp) site
#   (protospacer) in the human EMX1 locus that precedes an NGG trinucleotide,
#   the requisite protospacer-adjacent motif (PAM)".  CITED, NOT VENDORED --
#   the PMC copy is text-mining/fair-use only, so the NGG half of this grammar
#   is byte-pinned via CRISPOR and the paper is the human-readable source.
#
# AsCas12a: TTTV on the protospacer's 5' side, 23 nt guide.
#   Zetsche et al. 2015 (PMID 26422227, PMC4638220): "the PAM for FnCpf1 is
#   located upstream of the 5' end of the displaced strand of the protospacer
#   and has the sequence 5'-TTN".  CITED, NOT VENDORED -- see register CR-4.
#   The narrowing from TTTN to TTTV is a design convention, not a claim about
#   what the nuclease cuts, and the subject is held only to its own TTTV claim.
#
# LwaCas13a: no PAM, no second strand.
#   Abudayyeh et al. 2017 (PMID 28976959, PMC5706658): "motif analysis of the
#   depleted PFS sequences at varying thresholds revealed the expected 3' H
#   motif of LshCas13a, but no significant PFS motif for LwaCas13a".  CITED,
#   NOT VENDORED.  The subject claims no PFS and that claim agrees with the
#   literature, so CR-4's Cas13 clause tests self-consistency only.

CANONICAL_LENGTH = {"SpCas9": 20, "SpCas9-HF1": 20, "AsCas12a": 23, "LwaCas13a": 22}
PAM_MOTIF = {"SpCas9": "NGG", "SpCas9-HF1": "NGG", "AsCas12a": "TTTV", "LwaCas13a": "none"}
PAM_SIDE = {"SpCas9": 3, "SpCas9-HF1": 3, "AsCas12a": 5, "LwaCas13a": 0}
PAM_REGEX = {"NGG": r"[ACGTN]GG", "TTTV": r"TTT[ACG]"}


def rc(seq: str) -> str:
    """Biopython's reverse complement -- the oracle, not the subject's own."""
    return reverse_complement(seq)


def gc_percent(seq: str) -> float:
    """Biopython's GC fraction as a percentage, rounded the way the app rounds."""
    return round(gc_fraction(seq, ambiguous="remove") * 100, 1)


# ── independent both-strand candidate scans ──────────────────────────────────

def scan_cas9(seq: str, guide_length: int) -> list[tuple[str, int, str, str]]:
    """Every (strand, plus-strand start, protospacer, PAM) for an NGG nuclease.

    The plus-strand scan walks ``seq``; the minus-strand scan walks Biopython's
    reverse complement and maps each hit back by ``len - start - L``. Overlapping
    sites are kept -- a lookahead, not ``re.findall``.
    """
    n, L = len(seq), guide_length
    out = []
    pat = rf"(?=([ACGTN]{{{L}}})([ACGTN]GG))"
    for m in re.finditer(pat, seq):
        out.append(("+", m.start(), m.group(1), m.group(2)))
    for m in re.finditer(pat, rc(seq)):
        out.append(("-", n - m.start() - L, m.group(1), m.group(2)))
    return out


def scan_cas12a(seq: str, guide_length: int) -> list[tuple[str, int, str, str]]:
    """Every (strand, plus-strand start, protospacer, PAM) for a 5' TTTV PAM."""
    n, L = len(seq), guide_length
    out = []
    pat = rf"(?=(TTT[ACG])([ACGTN]{{{L}}}))"
    for m in re.finditer(pat, seq):
        out.append(("+", m.start() + 4, m.group(2), m.group(1)))
    for m in re.finditer(pat, rc(seq)):
        out.append(("-", n - (m.start() + 4) - L, m.group(2), m.group(1)))
    return out


def scan_cas13(seq: str, guide_length: int) -> list[tuple[str, int, str, str]]:
    """Every window of the given sense strand; no PAM, no second strand."""
    L = guide_length
    return [("+", i, seq[i:i + L], "none") for i in range(len(seq) - L + 1)]


SCAN = {
    "SpCas9": scan_cas9,
    "SpCas9-HF1": scan_cas9,
    "AsCas12a": scan_cas12a,
    "LwaCas13a": scan_cas13,
}


def scan(cas: str, seq: str, guide_length: int | None = None):
    return SCAN[cas](seq, guide_length or CANONICAL_LENGTH[cas])


def scan_circular(cas: str, seq: str, guide_length: int | None = None) -> int:
    """How many candidates the same grammar finds if the molecule is a circle.

    Every start on the circle is tried against the doubled sequence, so sites
    that run through the origin are counted exactly once each.
    """
    L = guide_length or CANONICAL_LENGTH[cas]
    n = len(seq)
    if cas == "LwaCas13a":
        return n
    span, pat = (L + 3, rf"[ACGTN]{{{L}}}[ACGTN]GG") if PAM_SIDE[cas] == 3 \
        else (L + 4, rf"TTT[ACG][ACGTN]{{{L}}}")
    total = 0
    for strand_seq in (seq, rc(seq)):
        doubled = strand_seq * 2
        total += sum(1 for i in range(n) if re.match(pat, doubled[i:i + span]))
    return total


# ── re-derivation by index: what the raw sequence says is at a reported site ──

def site_at(cas: str, seq: str, position: int, guide_length: int) -> tuple[str, str]:
    """(protospacer, PAM) read straight out of ``seq`` at a reported position.

    ``position`` is the leftmost plus-strand index of the protospacer footprint,
    for either strand -- the convention the subject reports. A minus-strand site
    is read by reverse-complementing the plus-strand slice, so the PAM comes from
    the opposite side of the footprint.
    """
    L, n = guide_length, len(seq)
    proto_fwd = seq[position:position + L]
    if cas == "LwaCas13a":
        return proto_fwd, "none"
    if PAM_SIDE[cas] == 3:                       # NGG, 3' of the protospacer
        plus = (proto_fwd, seq[position + L:position + L + 3])
        minus = (rc(proto_fwd), rc(seq[max(position - 3, 0):position]))
    else:                                        # TTTV, 5' of the protospacer
        plus = (proto_fwd, seq[max(position - 4, 0):position])
        minus = (rc(proto_fwd), rc(seq[position + L:min(position + L + 4, n)]))
    return plus, minus


def site_on(cas: str, seq: str, position: int, strand: str, guide_length: int):
    plus, minus = site_at(cas, seq, position, guide_length)
    if cas == "LwaCas13a":
        return plus, minus
    return plus if strand == "+" else minus


# ── Doench et al. 2014 Rule Set 1, re-typed from the pinned CRISPOR copy ─────
# gx/corpus/crispr/crisporEffScores-v3.1.py#doenchParams = [
# Positions are 0-based into the 30mer the model is defined on:
# gx/corpus/crispr/crisporEffScores-v3.1.py#Input is a 30mer: 4bp 5', 20bp guide, 3bp PAM, 3bp 5'

DOENCH_INTERCEPT = 0.59763615
DOENCH_GC_HIGH = -0.1665878
DOENCH_GC_LOW = -0.2026259
DOENCH_PARAMS = [
    (1, "G", -0.2753771), (2, "A", -0.3238875), (2, "C", 0.17212887), (3, "C", -0.1006662),
    (4, "C", -0.2018029), (4, "G", 0.24595663), (5, "A", 0.03644004), (5, "C", 0.09837684),
    (6, "C", -0.7411813), (6, "G", -0.3932644), (11, "A", -0.466099), (14, "A", 0.08537695),
    (14, "C", -0.013814), (15, "A", 0.27262051), (15, "C", -0.1190226), (15, "T", -0.2859442),
    (16, "A", 0.09745459), (16, "G", -0.1755462), (17, "C", -0.3457955), (17, "G", -0.6780964),
    (18, "A", 0.22508903), (18, "C", -0.5077941), (19, "G", -0.4173736), (19, "T", -0.054307),
    (20, "G", 0.37989937), (20, "T", -0.0907126), (21, "C", 0.05782332), (21, "T", -0.5305673),
    (22, "T", -0.8770074), (23, "C", -0.8762358), (23, "G", 0.27891626), (23, "T", -0.4031022),
    (24, "A", -0.0773007), (24, "C", 0.28793562), (24, "T", -0.2216372), (27, "G", -0.6890167),
    (27, "T", 0.11787758), (28, "C", -0.1604453), (29, "G", 0.38634258), (1, "GT", -0.6257787),
    (4, "GC", 0.30004332), (5, "AA", -0.8348362), (5, "TA", 0.76062777), (6, "GG", -0.4908167),
    (11, "GG", -1.5169074), (11, "TA", 0.7092612), (11, "TC", 0.49629861), (11, "TT", -0.5868739),
    (12, "GG", -0.3345637), (13, "GA", 0.76384993), (13, "GC", -0.5370252), (16, "TG", -0.7981461),
    (18, "GG", -0.6668087), (18, "TC", 0.35318325), (19, "CC", 0.74807209), (19, "TG", -0.3672668),
    (20, "AC", 0.56820913), (20, "CG", 0.32907207), (20, "GA", -0.8364568), (20, "GG", -0.7822076),
    (21, "TC", -1.029693), (22, "CG", 0.85619782), (22, "CT", -0.4632077), (23, "AA", -0.5794924),
    (23, "AG", 0.64907554), (24, "AG", -0.0773007), (24, "CG", 0.28793562), (24, "TG", -0.2216372),
    (26, "GT", 0.11787758), (28, "GG", -0.69774),
]


def doench_rule_set1(mer30: str) -> float:
    """Doench 2014 Rule Set 1 on a 30mer, in [0, 1].

    30mer layout, from the pinned definition: 4 nt of 5' context, the 20 nt
    protospacer, the 3 nt PAM, 3 nt of 3' context. The PAM and both flanks carry
    weights -- positions 24 and above in ``DOENCH_PARAMS`` are PAM/3'-context
    terms -- so the score is NOT a function of the protospacer alone.
    """
    if len(mer30) != 30:
        raise ValueError(f"Rule Set 1 is defined on a 30mer, got {len(mer30)}")
    score = DOENCH_INTERCEPT
    guide = mer30[4:24]
    gc = guide.count("G") + guide.count("C")
    score += abs(10 - gc) * (DOENCH_GC_LOW if gc <= 10 else DOENCH_GC_HIGH)
    for pos, model_seq, weight in DOENCH_PARAMS:
        if mer30[pos:pos + len(model_seq)] == model_seq:
            score += weight
    return 1.0 / (1.0 + math.exp(-score))


def mer30_for(seq: str, position: int, strand: str, guide_length: int = 20) -> str | None:
    """The Rule Set 1 input window around a reported SpCas9 site, or None.

    None when the site sits too close to an end of the molecule for the 4 nt of
    5' context or the 3 nt of 3' context to exist -- Rule Set 1 is undefined
    there, and guessing the missing bases would be inventing an expected value.
    """
    if strand == "-":
        seq, position = rc(seq), len(seq) - position - guide_length
    start, end = position - 4, position + guide_length + 6
    if start < 0 or end > len(seq):
        return None
    return seq[start:end]


# ── Hsu et al. 2013 / MIT, re-typed from the pinned CRISPOR copy ─────────────
# gx/corpus/crispr/crispor-v3.1.py#hitScoreM = [0,0,0.014,...]
# gx/corpus/crispr/crispor-v3.1.py#def calcHitScore(string1,string2):
# gx/corpus/crispr/crispor-v3.1.py#def calcMitGuideScore(hitSum):

HIT_SCORE_M = [
    0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583,
]


def mit_hit_score(guide: str, candidate: str) -> float:
    """MIT single-hit score in [0, 100], over the ACTUAL mismatch positions.

    Transcribed from ``calcHitScore``. The whole point of the published score is
    that it is position-dependent: index 19 (protospacer position 20, in the
    PAM-proximal seed) carries weight 0.583, index 0 (the PAM-distal 5' base)
    carries 0.
    """
    assert len(guide) == len(candidate) == 20, "transcribed for the 20mer case"
    dists, mm_count, last = [], 0, None
    score1 = 1.0
    for pos in range(20):
        if guide[pos] != candidate[pos]:
            mm_count += 1
            if last is not None:
                dists.append(pos - last)
            score1 *= 1 - HIT_SCORE_M[pos]
            last = pos
    score2 = 1.0 if mm_count < 2 else 1.0 / (((19 - sum(dists) / len(dists)) / 19.0) * 4 + 1)
    score3 = 1.0 if mm_count == 0 else 1.0 / (mm_count ** 2)
    return score1 * score2 * score3 * 100


def mit_guide_score(hit_sum: float) -> float:
    """``calcMitGuideScore`` without its final int-rounding: 100/(100+hitSum)*100."""
    return 100.0 / (100.0 + hit_sum) * 100


# ── the re-typed tables must still match the pinned bytes ────────────────────

def _assert_against_pinned_bytes() -> None:
    eff = (CORPUS / "crispr" / "crisporEffScores-v3.1.py").read_text()
    main = (CORPUS / "crispr" / "crispor-v3.1.py").read_text()

    pinned_params = eval("[" + eff.split("doenchParams = [")[1].split("]")[0] + "]")
    assert pinned_params == DOENCH_PARAMS, "DOENCH_PARAMS drifted from the pinned bytes"
    for literal in (f"intercept =  {DOENCH_INTERCEPT}",
                    f"gcHigh    = {DOENCH_GC_HIGH}",
                    f"gcLow     = {DOENCH_GC_LOW}",
                    "Input is a 30mer: 4bp 5', 20bp guide, 3bp PAM, 3bp 5'"):
        assert literal in eff, f"anchor missing from pinned bytes: {literal!r}"

    pinned_m = eval(main.split("hitScoreM = ")[1].split("]")[0] + "]")
    assert pinned_m == HIT_SCORE_M, "HIT_SCORE_M drifted from the pinned bytes"
    for literal in ("def calcHitScore(string1,string2):",
                    "def calcMitGuideScore(hitSum):",
                    "('NGG','NGG - Streptococcus Pyogenes')"):
        assert literal in main, f"anchor missing from pinned bytes: {literal!r}"


_assert_against_pinned_bytes()


# ── small helpers the tests share ────────────────────────────────────────────

def spearman(xs: list[float], ys: list[float]) -> float:
    """Rank correlation, ties broken by position -- enough to rank 300+ guides."""
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        out = [0] * len(v)
        for k, i in enumerate(order):
            out[i] = k
        return out
    a, b = ranks(xs), ranks(ys)
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return num / den if den else 0.0


def overlaps(a_pos: int, b_pos: int, length: int) -> bool:
    return abs(a_pos - b_pos) < length
