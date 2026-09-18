"""PCR primer design and restriction enzyme analysis."""

import re
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from Bio import Restriction
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp

router = APIRouter(prefix="/primers", tags=["primers"])


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


class Primer(BaseModel):
    sequence: str
    position: int
    length: int
    tm: float
    gc_content: float
    direction: str  # forward | reverse


class PrimerPair(BaseModel):
    forward: Primer
    reverse: Primer
    product_size: int
    penalty: float


class RestrictionSite(BaseModel):
    enzyme: str
    cut_pattern: str
    positions: list[int]
    count: int


#: Conditions the reported Tm is computed under. Spelled out because a melting
#: temperature is not a property of a sequence alone: Tm_NN takes a
#: nearest-neighbour table, both strand concentrations, four salt species and a
#: choice of salt correction, and "Tm" without them does not identify a number.
#: These are Biopython's documented defaults, which is what the #50 audit compared
#: against. Changing any of them changes every Tm the app reports, so change them
#: here, visibly, and update backend/tests/test_oracle_tm.py in the same commit.
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


def _tm_nearest_neighbor(seq: str) -> float:
    """Nearest-neighbour melting temperature, as the name says.

    Delegates to Biopython's `Tm_NN` under the pinned `TM_CONDITIONS`. An
    independent audit of the previous Wallace-rule implementation found 97 of 120
    primers from real templates more than 3 degC out, worst case +9.0 degC, and
    systematically worse on GC-rich sequence — which is where an annealing
    temperature set from this number does real damage (#50).
    """
    seq = seq.upper()
    if len(seq) < 2:
        # Tm_NN needs at least one neighbour pair; nothing thermodynamic to say.
        return _tm_wallace(seq)
    return round(MeltingTemp.Tm_NN(seq, **TM_CONDITIONS), 2)


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

    template = req.template.upper().replace(" ", "").replace("\n", "")
    if len(template) < req.product_min + req.primer_len_min * 2:
        raise HTTPException(400, "Template too short for requested product size")

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
            tm = _tm_nearest_neighbor(seq)
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
            tm = _tm_nearest_neighbor(rc_seq)
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

            counter += 1
            entry = (-penalty, counter, fwd, rev, rev.position + rev.length - f_pos)
            if len(best) < req.max_pairs:
                heapq.heappush(best, entry)
                if len(best) == req.max_pairs:
                    limit = min(5.0, -best[0][0])
            elif penalty < -best[0][0]:
                heapq.heapreplace(best, entry)
                limit = min(5.0, -best[0][0])

    pairs = [
        PrimerPair(forward=f, reverse=r, product_size=prod, penalty=round(-neg, 3))
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
    template = template.upper().replace(" ", "").replace("\n", "")  # noqa: F841 (reassigned from req)
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
    template = req.template.upper().replace(" ", "").replace("\n", "")
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
