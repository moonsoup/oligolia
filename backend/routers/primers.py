"""PCR primer design and restriction enzyme analysis."""

import math
import re
from typing import Literal
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from Bio import Restriction
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp

from backend.sequence_text import normalize_template

router = APIRouter(prefix="/primers", tags=["primers"])


#: Conditions the reported Tm is computed under. Spelled out because a melting
#: temperature is not a property of a sequence alone: Tm_NN takes a
#: nearest-neighbour table, both strand concentrations, four salt species and a
#: choice of salt correction, and "Tm" without them does not identify a number.
#: These are Biopython's documented defaults, which is what the #50 audit compared
#: against. Changing any of them changes every Tm the app reports, so change them
#: here, visibly, and update backend/tests/test_oracle_tm.py in the same commit.
#:
#: Since #81 a caller may ASK for a different buffer per request, but these remain
#: the defaults, unchanged: whether they are the right defaults for a PCR tool is
#: #79's open scientific question and is deliberately not settled here.
TM_CONDITIONS = {
    "nn_table": MeltingTemp.DNA_NN3,  # Allawi & SantaLucia 1997 ("unified")
    "saltcorr": 5,                    # Owczarzy et al. 2004
    "dnac1": 25,                      # nM, primer strand
    "dnac2": 25,                      # nM, template strand
    "Na": 50,                         # mM
    "K": 0,
    "Tris": 0,
    "Mg": 0,
    "dNTPs": 0,
}


class PrimerDesignRequest(BaseModel):
    template: str
    product_min: int = Field(default=100, ge=50)
    product_max: int = Field(default=1000, le=10000)
    primer_len_min: int = Field(default=18, ge=15)
    primer_len_max: int = Field(default=24, le=35)
    tm_min: float = 55.0
    tm_max: float = 65.0
    gc_min: float = 40.0
    gc_max: float = 70.0
    max_pairs: int = Field(default=5, ge=1, le=20)

    # --- how to compute Tm (#81) ---------------------------------------------
    # "nn" is nearest-neighbour thermodynamics via Tm_NN (what the app reports and
    # what #50 established it must report). "wallace" is the 1989 length-and-GC-count
    # rule, kept since #50 and now actually selectable — #50 called it "an option"
    # while nothing could reach it, which is what #79.2 objected to.
    tm_method: Literal["nn", "wallace"] = "nn"

    # The buffer the Tm is computed in. Defaults are EXACTLY TM_CONDITIONS, so a
    # request that sets none of these is byte-for-byte the pre-#81 request.
    # ge=0 rather than gt=0: zero is a legitimate concentration (the defaults have
    # four of them), but a negative one is not a buffer, so it is a 422.
    na_mm: float = Field(default=TM_CONDITIONS["Na"], ge=0)
    k_mm: float = Field(default=TM_CONDITIONS["K"], ge=0)
    tris_mm: float = Field(default=TM_CONDITIONS["Tris"], ge=0)
    mg_mm: float = Field(default=TM_CONDITIONS["Mg"], ge=0)
    dntp_mm: float = Field(default=TM_CONDITIONS["dNTPs"], ge=0)
    primer_nm: float = Field(default=TM_CONDITIONS["dnac1"], ge=0)
    template_nm: float = Field(default=TM_CONDITIONS["dnac2"], ge=0)


class Primer(BaseModel):
    sequence: str
    position: int
    length: int
    tm: float
    gc_content: float
    direction: str  # forward | reverse


class TmConditions(BaseModel):
    """The buffer a reported Tm was computed in (#81).

    Carried on every pair because a Tm without its conditions does not identify a
    number — the same primer reads 53.7 degC in the default bench buffer and 64.0
    degC in a PCR-like one (#79). A caller that never asks for a buffer still gets
    this, so a reported Tm is never again silently conditioned.
    """
    na_mm: float
    k_mm: float
    tris_mm: float
    mg_mm: float
    dntp_mm: float
    primer_nm: float
    template_nm: float
    # Only meaningful for the nearest-neighbour method; None under the Wallace
    # rule, which is a function of length and GC count and sees no buffer at all.
    nn_table: str | None = None
    salt_correction: int | None = None


class PrimerPair(BaseModel):
    forward: Primer
    reverse: Primer
    product_size: int
    penalty: float
    # Additive (#81): the response stays a JSON list of pairs, so every existing
    # caller — the GUI, workflow/engine.py, the QA corpus — is unaffected.
    tm_method: str = "nn"
    tm_conditions: TmConditions | None = None


class RestrictionSite(BaseModel):
    enzyme: str
    cut_pattern: str
    positions: list[int]
    count: int


def _tm_wallace(seq: str) -> float:
    """The original approximation, kept available and no longer misnamed.

    This is what `_tm_nearest_neighbor` used to compute: the 1989 Wallace rule
    (2*AT + 4*GC) below 14 nt, and the Marmur/Chester GC-percentage formula at or
    above it. Both are functions of length and GC *count* only, so they cannot see
    adjacent-base stacking — two primers of equal length and GC count get the same
    answer whatever order the bases are in.

    Retained because it is cheap and is a reasonable pre-filter, and because
    removing a function is not this project's call to make unilaterally. It is not
    what the app reports (#50).
    """
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    if len(seq) < 14:
        return 2 * at + 4 * gc
    return 64.9 + 41 * (gc - 16.4) / (at + gc)


def _tm_nearest_neighbor(seq: str, conditions: dict | None = None) -> float:
    """Nearest-neighbour melting temperature, as the name says.

    Delegates to Biopython's `Tm_NN`. An independent audit of the previous
    Wallace-rule implementation found 97 of 120 primers from real templates more
    than 3 degC out, worst case +9.0 degC, and systematically worse on GC-rich
    sequence — which is where an annealing temperature set from this number does
    real damage (#50).

    `conditions` is a Tm_NN keyword dict; omitting it uses the pinned
    `TM_CONDITIONS`, so every existing caller keeps its exact previous behaviour
    (#81 adds the argument rather than changing the default).
    """
    seq = seq.upper()
    if len(seq) < 2:
        # Tm_NN needs at least one neighbour pair; nothing thermodynamic to say.
        return _tm_wallace(seq)
    return round(MeltingTemp.Tm_NN(seq, **(TM_CONDITIONS if conditions is None else conditions)), 2)


def _effective_monovalent_mm(conditions: dict) -> float:
    """The total monovalent-equivalent ion concentration Tm_NN's saltcorr will see.

    Mirrors `Bio.SeqUtils.MeltingTemp.salt_correction`: monovalents add directly,
    Tris counts half (only half of it is the cation), and free Mg2+ — magnesium
    not chelated by dNTPs — enters as a von Ahsen (2001) sodium equivalent. If the
    whole thing is zero, corrections 1-6 raise rather than return a number.
    """
    mon = conditions["Na"] + conditions["K"] + conditions["Tris"] / 2.0
    mg, dntps = conditions["Mg"], conditions["dNTPs"]
    if sum((conditions["K"], mg, conditions["Tris"], dntps)) > 0 and dntps < mg:
        mon += 120 * math.sqrt(mg - dntps)
    return mon


def _validate_tm_conditions(conditions: dict) -> None:
    """Refuse buffers Tm_NN cannot evaluate, with a 4xx rather than a 500 (#81).

    Both of these are `ValueError` out of Biopython — one raised deliberately, one
    a bare `math domain error` — and both surface as an unexplained 500 if they
    reach the endpoint. They are properties of the requested buffer, so they are
    the caller's to fix and belong in the 4xx range.
    """
    if _effective_monovalent_mm(conditions) <= 0:
        raise HTTPException(400, (
            "No monovalent salt and no free Mg²⁺: the salt correction takes the "
            "logarithm of the total ion concentration, which is zero here. Give a "
            "non-zero na_mm, k_mm or tris_mm, or mg_mm greater than dntp_mm."
        ))
    if conditions["dnac1"] - conditions["dnac2"] / 2.0 <= 0:
        raise HTTPException(400, (
            f"primer_nm ({conditions['dnac1']}) must be greater than half of "
            f"template_nm ({conditions['dnac2']}): the nearest-neighbour "
            "concentration term takes the logarithm of primer − template/2."
        ))


def _tm_conditions_for(req: "PrimerDesignRequest") -> dict:
    """Tm_NN keywords for this request — the pinned table and correction, its buffer."""
    return {
        "nn_table": TM_CONDITIONS["nn_table"],
        "saltcorr": TM_CONDITIONS["saltcorr"],
        "dnac1": req.primer_nm,
        "dnac2": req.template_nm,
        "Na": req.na_mm,
        "K": req.k_mm,
        "Tris": req.tris_mm,
        "Mg": req.mg_mm,
        "dNTPs": req.dntp_mm,
    }


def _nn_table_name(table: dict) -> str:
    """Biopython's nearest-neighbour tables are plain dicts, so ask it their name.

    Reported rather than assumed: "Tm 58.8 degC" means nothing without saying which
    of the four literature parameter sets produced it.
    """
    for name in ("DNA_NN1", "DNA_NN2", "DNA_NN3", "DNA_NN4", "RNA_NN1",
                 "RNA_NN2", "RNA_NN3", "R_DNA_NN1"):
        if getattr(MeltingTemp, name, None) is table:
            return name
    return "unknown"


def _reported_conditions(conditions: dict, method: str) -> TmConditions:
    """The buffer, as it goes back to the caller on every pair."""
    return TmConditions(
        na_mm=conditions["Na"],
        k_mm=conditions["K"],
        tris_mm=conditions["Tris"],
        mg_mm=conditions["Mg"],
        dntp_mm=conditions["dNTPs"],
        primer_nm=conditions["dnac1"],
        template_nm=conditions["dnac2"],
        # The Wallace rule cannot see a nearest-neighbour table or a salt
        # correction, so claiming one would be a false report.
        nn_table=_nn_table_name(conditions["nn_table"]) if method == "nn" else None,
        salt_correction=conditions["saltcorr"] if method == "nn" else None,
    )


#: Human-readable names for the pinned table and salt correction, for the one
#: place a user sees them: the GUI's Tm column caption (#81, #79.4).
TM_TABLE_NAME = "DNA_NN3 (Allawi & SantaLucia 1997)"
TM_SALTCORR_NAME = "Owczarzy et al. 2004"


def describe_tm_conditions(conditions: dict | None = None, method: str = "nn") -> str:
    """One line naming the Tm method and the buffer it assumes.

    The primers panel shows a Tm column; before #81 it said nothing about what
    that number assumes, which is exactly the complaint in #79.4. Derived from
    `TM_CONDITIONS` rather than retyped, so the caption cannot drift away from the
    number beside it.
    """
    c = TM_CONDITIONS if conditions is None else conditions
    buffer = (
        f"Na⁺ {c['Na']:g} mM, K⁺ {c['K']:g} mM, Tris {c['Tris']:g} mM, "
        f"Mg²⁺ {c['Mg']:g} mM, dNTPs {c['dNTPs']:g} mM, "
        f"primer {c['dnac1']:g} nM, template {c['dnac2']:g} nM"
    )
    if method == "wallace":
        return ("Tm: Wallace rule (length and GC count only) — ignores the buffer; "
                f"the reaction is assumed to be {buffer}.")
    return (f"Tm: nearest-neighbour, {TM_TABLE_NAME}, salt correction "
            f"{TM_SALTCORR_NAME}, at {buffer}.")


def describe_pair_tm(pair: "PrimerPair") -> str:
    """The same line, for a pair that has already been designed.

    Lets a table caption report what its rows were ACTUALLY computed under rather
    than what the defaults happen to be.
    """
    c = pair.tm_conditions
    if c is None:
        return describe_tm_conditions(method=pair.tm_method)
    return describe_tm_conditions({
        "Na": c.na_mm, "K": c.k_mm, "Tris": c.tris_mm,
        "Mg": c.mg_mm, "dNTPs": c.dntp_mm,
        "dnac1": c.primer_nm, "dnac2": c.template_nm,
    }, method=pair.tm_method)


def _gc(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq) * 100 if seq else 0.0


#: The shortest loop a primer can actually fold around. Below 3 nt the backbone
#: cannot turn, so two complementary arms sitting closer than this are not a
#: hairpin however well they pair.
MIN_HAIRPIN_LOOP = 3

#: Stem length that counts as a hairpin. At 4 bp even ordinary sequence pairs up
#: somewhere — a 4-mer has only 256 possibilities — so a 4-bp stem with a loop is
#: not worth rejecting a primer over. Real tools score hairpins by free energy;
#: this is a length heuristic and is deliberately set where it does not reject
#: most of the template.
MIN_HAIRPIN_STEM = 5

#: At most this many of the last three bases may be G or C. A strong 3' end
#: promotes mispriming.
MAX_3PRIME_GC = 2

#: Homopolymer runs to reject outright.
HOMOPOLYMERS = ("AAAA", "TTTT", "GGGG", "CCCC")


def _has_hairpin(seq: str, min_stem: int = MIN_HAIRPIN_STEM, min_loop: int = MIN_HAIRPIN_LOOP) -> bool:
    """Does the primer fold back on itself?

    A hairpin needs two things: a stem of complementary arms, AND a loop between
    them long enough for the backbone to turn. The previous version checked only
    whether any `min_stem`-mer appeared anywhere in the whole reverse complement,
    with no loop requirement and no positional relationship at all — so it flagged
    56% of random 20-mers.

    `GAATTCAAAAAAAAAAAAAA` was its headline false positive: GAATTC is
    self-complementary, but its two halves are ADJACENT, so the loop length is
    zero and it cannot fold (#62).
    """
    seq = seq.upper()
    comp = str.maketrans("ACGT", "TGCA")
    n = len(seq)

    for i in range(n - 2 * min_stem - min_loop + 1):
        arm = seq[i:i + min_stem]
        # The downstream arm has to be this arm's reverse complement, and it has
        # to start at least min_loop bases after the upstream arm ends.
        target = arm.translate(comp)[::-1]
        if target in seq[i + min_stem + min_loop:]:
            return True
    return False


def _acceptable_primer(seq: str) -> bool:
    """The sequence-only rules a primer must pass, whichever orientation it is.

    Forward and reverse candidates were filtered by DIFFERENT rules: reverse
    candidates skipped the 3'-GC-clamp and homopolymer checks entirely, so 14 of
    20 returned pairs had reverse primers that the forward filter would have
    rejected — `GCTGCCCCTACGGATCGCA` among them (#62).

    Tm and GC window stay at the call site because they come from the request.
    """
    seq = seq.upper()
    if _has_hairpin(seq):
        return False
    if seq[-3:].count("G") + seq[-3:].count("C") > MAX_3PRIME_GC:
        return False
    return not any(run in seq for run in HOMOPOLYMERS)


#: The longest 3'-end complementarity a pair may have. Three or fewer bases pair
#: by chance often enough that rejecting on them would throw away most usable
#: pairs -- the mistake #62's hairpin filter made.
MAX_3PRIME_DIMER = 3


def three_prime_dimer_length(fwd: str, rev: str) -> int:
    """How many bases of the two primers' 3' ends are complementary.

    A 3' primer-dimer is a property of the PAIR, not of either primer, so no
    amount of per-primer filtering finds it (#52). It is one of the commonest
    reasons a PCR consumes itself without amplifying the target: the two primers
    anneal to each other by their 3' ends and extend each other instead.

    Only the 3' ends matter. Complementarity at the 5' ends cannot prime, because
    it is the 3' end that a polymerase extends -- so this compares the tail of one
    primer against the reverse complement of the tail of the other, and returns
    the longest such overlap found.
    """
    fwd, rev = fwd.upper(), rev.upper()
    comp = str.maketrans("ACGT", "TGCA")
    longest = 0

    for n in range(1, min(len(fwd), len(rev)) + 1):
        tail = fwd[-n:]
        # If rev's 3' tail is the reverse complement of fwd's 3' tail, they anneal
        # 3'-to-3' and both can extend.
        if rev[-n:] == tail.translate(comp)[::-1]:
            longest = n
    return longest


def forms_3prime_dimer(fwd: str, rev: str, limit: int = MAX_3PRIME_DIMER) -> bool:
    """Is the pair's 3'-end complementarity long enough to matter?"""
    return three_prime_dimer_length(fwd, rev) > limit


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


@router.post("/design", response_model=list[PrimerPair])
def design_primers(req: PrimerDesignRequest) -> list[PrimerPair]:
    if req.tm_min > req.tm_max:
        raise HTTPException(400, f"tm_min ({req.tm_min}°C) must be ≤ tm_max ({req.tm_max}°C)")
    if req.product_min > req.product_max:
        raise HTTPException(400, f"product_min ({req.product_min}) must be ≤ product_max ({req.product_max})")
    if req.primer_len_min > req.primer_len_max:
        raise HTTPException(400, f"primer_len_min ({req.primer_len_min}) must be ≤ primer_len_max ({req.primer_len_max})")
    if req.gc_min > req.gc_max:
        raise HTTPException(400, f"gc_min ({req.gc_min}) must be ≤ gc_max ({req.gc_max})")

    template = normalize_template(req.template)
    if len(template) < req.product_min + req.primer_len_min * 2:
        raise HTTPException(400, "Template too short for requested product size")

    # The requested buffer, checked before any Tm is computed (#81). Validated
    # whichever method was asked for, so a refusal is a property of the buffer and
    # not of the method: the conditions go back on every pair either way, and a
    # buffer we could not report a Tm under is not one to echo back as if it were
    # fine.
    conditions = _tm_conditions_for(req)
    _validate_tm_conditions(conditions)
    reported_conditions = _reported_conditions(conditions, req.tm_method)

    if req.tm_method == "wallace":
        def _tm(seq: str) -> float:
            return _tm_wallace(seq)
    else:
        def _tm(seq: str) -> float:
            return _tm_nearest_neighbor(seq, conditions)

    fwd_candidates: list[Primer] = []
    rev_candidates: list[Primer] = []

    # Order matters for cost: GC and the sequence-only rules are cheap string
    # work, while _tm_nearest_neighbor is a nearest-neighbour calculation over
    # every dinucleotide. Computing Tm for all ~21k candidates of a 3 kb template
    # and then discarding most of them is what made this slow once #50 replaced
    # the old arithmetic Wallace formula (#62.3).
    for length in range(req.primer_len_min, req.primer_len_max + 1):
        for pos in range(0, len(template) - length + 1):
            seq = template[pos:pos + length]
            gc = _gc(seq)
            if not (req.gc_min <= gc <= req.gc_max):
                continue
            if not _acceptable_primer(seq):
                continue
            tm = _tm(seq)
            if not (req.tm_min <= tm <= req.tm_max):
                continue
            fwd_candidates.append(Primer(
                sequence=seq, position=pos, length=length,
                tm=round(tm, 1), gc_content=round(gc, 1), direction="forward",
            ))

        for pos in range(length, len(template) + 1):
            rc_seq = _reverse_complement(template[pos - length:pos])
            gc = _gc(rc_seq)
            if not (req.gc_min <= gc <= req.gc_max):
                continue
            # The same rules as forward. Reverse candidates used to skip the
            # 3'-GC-clamp and homopolymer checks entirely (#62.2).
            if not _acceptable_primer(rc_seq):
                continue
            tm = _tm(rc_seq)
            if not (req.tm_min <= tm <= req.tm_max):
                continue
            rev_candidates.append(Primer(
                sequence=rc_seq, position=pos - length, length=length,
                tm=round(tm, 1), gc_content=round(gc, 1), direction="reverse",
            ))

    # Index the reverse candidates by their 3' end so each forward primer only
    # visits the ones that can give a legal product. The old loop was O(F x R):
    # 14.5k x 14.5k on a 3 kb template, and it built a pydantic PrimerPair for
    # every one of the ~14.5M legal combinations before sorting and throwing all
    # but `max_pairs` away (#62.3).
    #
    # Two changes: bisect to the legal window, and keep only the best max_pairs
    # in a bounded heap of plain tuples, constructing PrimerPair objects once at
    # the end. Same pairs, same order, without materialising the cross product.
    import bisect
    import heapq

    rev_by_end = sorted(rev_candidates, key=lambda r: r.position + r.length)
    rev_ends = [r.position + r.length for r in rev_by_end]

    # A max-heap keyed on penalty (negated), so the worst kept pair is at the top
    # and can be evicted in O(log n).
    best: list[tuple[float, int, Primer, Primer, int]] = []
    counter = 0

    # `penalty = tm_diff + 0.1 * gc_diff`, and gc_diff >= 0, so penalty >= tm_diff
    # ALWAYS. That makes tm_diff an admissible lower bound: once the heap is full,
    # any pair whose tm_diff already exceeds the worst kept penalty cannot beat it,
    # and can be rejected without computing the GC term or rounding anything.
    #
    # This is an exact prune, not a heuristic — the same pairs come out. It matters
    # because profiling showed the inner loop dominating everything else:
    # 24.7M round() and 60M abs() calls, against 1.1s total inside Tm_NN.
    limit = 5.0

    for fwd in fwd_candidates:
        f_pos, f_tm, f_gc = fwd.position, fwd.tm, fwd.gc_content
        lo = bisect.bisect_left(rev_ends, f_pos + req.product_min)
        hi = bisect.bisect_right(rev_ends, f_pos + req.product_max)

        for rev in rev_by_end[lo:hi]:
            tm_diff = f_tm - rev.tm
            if tm_diff < 0.0:
                tm_diff = -tm_diff
            if tm_diff > limit:
                continue
            if rev.position <= f_pos:
                continue

            gc_diff = f_gc - rev.gc_content
            if gc_diff < 0.0:
                gc_diff = -gc_diff
            penalty = tm_diff + gc_diff * 0.1

            full = len(best) >= req.max_pairs
            if full and penalty >= -best[0][0]:
                continue

            # A pair property, so it cannot be checked when the candidates are
            # built (#52). Deliberately LAST: it is string work, and by this point
            # the pair is both viable and good enough to be kept, so it runs on a
            # few hundred pairs rather than on millions.
            if forms_3prime_dimer(fwd.sequence, rev.sequence):
                continue

            counter += 1
            entry = (-penalty, counter, fwd, rev, rev.position + rev.length - f_pos)
            if not full:
                heapq.heappush(best, entry)
                if len(best) == req.max_pairs:
                    limit = min(5.0, -best[0][0])
            else:
                heapq.heapreplace(best, entry)
                limit = min(5.0, -best[0][0])

    pairs = [
        PrimerPair(forward=f, reverse=r, product_size=prod, penalty=round(-neg, 3),
                   tm_method=req.tm_method, tm_conditions=reported_conditions)
        for neg, _n, f, r, prod in sorted(best, key=lambda e: (-e[0], e[1]))
    ]
    return pairs


# Curated enzyme panel. Recognition sequences AND cut geometry (overhangs,
# real cut positions) are sourced from Bio.Restriction — the authoritative,
# vendored dataset — rather than hand-maintained (issue #34).
_CURATED_ENZYMES = [
    "EcoRI", "BamHI", "HindIII", "SalI", "XbaI", "SmaI", "KpnI", "SacI",
    "NotI", "XhoI", "NcoI", "NdeI", "ClaI", "SphI", "PstI", "PvuII",
    "AvaI", "AvaII", "EcoRV", "MluI", "NheI", "AgeI", "BglII", "MfeI",
]
_ENZYMES = {name: getattr(Restriction, name) for name in _CURATED_ENZYMES}
# name -> IUPAC recognition sequence (used for site mapping in restriction_sites)
RESTRICTION_ENZYMES = {name: enz.site for name, enz in _ENZYMES.items()}

IUPAC = {"R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
         "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
         "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}


def _pattern(recog: str) -> str:
    return "".join(IUPAC.get(c, c) for c in recog.upper())


def _find_sites(template: str, recog: str, is_circular: bool) -> list[int]:
    """0-indexed start positions of ``recog`` in ``template``.

    When ``is_circular`` is set, also finds recognition sites that span the
    origin junction (last k-1 bases + first k-1 bases, k = recognition length),
    reporting them at their real start index near the end of the template.
    """
    pattern = _pattern(recog)
    positions = [m.start() for m in re.finditer(f"(?={pattern})", template)]
    k, n = len(recog), len(template)
    if is_circular and k > 1 and n >= k:
        wrap = template[-(k - 1):] + template[:k - 1]
        for m in re.finditer(f"(?={pattern})", wrap):
            off = m.start()
            if off <= k - 2:  # starts in the tail => crosses the origin
                positions.append(n - (k - 1) + off)
    return sorted(set(positions))


class RestrictionRequest(BaseModel):
    template: str
    enzymes: list[str] | None = None
    is_circular: bool = False


@router.post("/restriction_sites", response_model=list[RestrictionSite])
def restriction_sites(req: RestrictionRequest) -> list[RestrictionSite]:
    template, enzymes = req.template, req.enzymes
    """Find restriction enzyme cut sites in a sequence."""
    # Normalised before any index is computed, so the positions reported here
    # are in the same coordinate space as /primers/digest's cuts (#86).
    template = normalize_template(template)  # noqa: F841 (reassigned from req)
    target_enzymes = {k: v for k, v in RESTRICTION_ENZYMES.items()
                      if not enzymes or k in enzymes}
    results = []
    for name, recog in target_enzymes.items():
        positions = _find_sites(template, recog, req.is_circular)
        if positions:
            results.append(RestrictionSite(
                enzyme=name,
                cut_pattern=recog,
                positions=positions,
                count=len(positions),
            ))
    results.sort(key=lambda r: (-r.count, r.enzyme))
    return results


class DigestRequest(BaseModel):
    template: str
    enzymes: list[str]
    is_circular: bool = False


class DigestFragment(BaseModel):
    start: int
    end: int
    length: int
    sequence: str
    # Single-stranded overhang produced by the cut at each end. Type is
    # "5'", "3'", "blunt" (blunt cutter), or "none" (a free linear terminus).
    # Overhang bases are given on the top strand; for the palindromic-site
    # enzymes in this panel that equals the complementary end's overhang.
    left_overhang: str = ""
    left_overhang_type: str = "none"
    right_overhang: str = ""
    right_overhang_type: str = "none"


def _overhang_at(template: str, cut: int, ovhg: int, is_circular: bool) -> tuple[str, str]:
    """Return (overhang_bases, type) for a top-strand cut at 0-based ``cut``.

    ``ovhg`` follows Bio.Restriction's sign convention: negative = 5' overhang,
    positive = 3' overhang, 0 = blunt.
    """
    if ovhg == 0:
        return "", "blunt"
    if ovhg < 0:  # 5' overhang spans [cut, cut - ovhg)
        idxs = range(cut, cut - ovhg)
        oh_type = "5'"
    else:  # 3' overhang spans [cut - ovhg, cut)
        idxs = range(cut - ovhg, cut)
        oh_type = "3'"
    n = len(template)
    if is_circular:
        bases = "".join(template[i % n] for i in idxs)
    else:
        bases = "".join(template[i] for i in idxs if 0 <= i < n)
    return bases, oh_type


class DigestResult(BaseModel):
    enzymes: list[str]
    cut_positions: list[int]
    fragments: list[DigestFragment]
    template_length: int


@router.post("/digest", response_model=DigestResult)
def digest(req: DigestRequest) -> DigestResult:
    """Simulate restriction digest — cut template at all sites for the given enzymes."""
    # Normalised the way Bio.Restriction normalises internally, so the string
    # sliced below is the string its cut positions index into (#86).
    template = normalize_template(req.template)
    if not template:
        raise HTTPException(400, "Template sequence is required")
    unknown = [e for e in req.enzymes if e not in RESTRICTION_ENZYMES]
    if unknown:
        raise HTTPException(400, f"Unknown enzymes: {unknown}. Supported: {sorted(RESTRICTION_ENZYMES)}")

    n = len(template)
    # Real cut geometry from Bio.Restriction: search() gives 1-based top-strand
    # cut positions; 0-based cut index (where the downstream fragment starts) is
    # position - 1. Overhang for each cut is derived from the enzyme's ovhg.
    bio_seq = Seq(template)
    linear = not req.is_circular
    cut_overhangs: dict[int, tuple[str, str]] = {}
    for enzyme in req.enzymes:
        enz = _ENZYMES[enzyme]
        for pos in enz.search(bio_seq, linear=linear):
            cut = (pos - 1) % n if req.is_circular else (pos - 1)
            cut_overhangs[cut] = _overhang_at(template, cut, enz.ovhg, req.is_circular)

    cut_positions = sorted(cut_overhangs)

    def _ends(start: int, end: int) -> dict:
        lo, lt = cut_overhangs.get(start, ("", "none"))
        ro, rt = cut_overhangs.get(end, ("", "none"))
        return {"left_overhang": lo, "left_overhang_type": lt,
                "right_overhang": ro, "right_overhang_type": rt}

    fragments: list[DigestFragment] = []
    if req.is_circular:
        # A circular molecule with N cuts yields N fragments (no free ends);
        # the fragment from the last cut wraps the origin back to the first.
        if not cut_positions:
            fragments.append(DigestFragment(start=0, end=n, length=n, sequence=template))
        else:
            count = len(cut_positions)
            for i in range(count):
                start = cut_positions[i]
                end = cut_positions[(i + 1) % count]
                if i < count - 1:
                    seq = template[start:end]
                    length = end - start
                else:  # origin-spanning fragment
                    seq = template[start:] + template[:end]
                    length = (n - start) + end
                fragments.append(DigestFragment(
                    start=start, end=end, length=length, sequence=seq, **_ends(start, end),
                ))
    else:
        # Linear: N cuts yield N+1 fragments, including the two end pieces.
        boundaries = [0] + cut_positions + [n]
        for i in range(len(boundaries) - 1):
            start, end = boundaries[i], boundaries[i + 1]
            seq = template[start:end]
            if seq:
                fragments.append(DigestFragment(
                    start=start, end=end, length=len(seq), sequence=seq, **_ends(start, end),
                ))

    fragments.sort(key=lambda f: -f.length)
    return DigestResult(
        enzymes=req.enzymes,
        cut_positions=cut_positions,
        fragments=fragments,
        template_length=len(template),
    )
