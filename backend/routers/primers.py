"""PCR primer design and restriction enzyme analysis."""

import re
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

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


def _tm_nearest_neighbor(seq: str) -> float:
    """Wallace rule Tm approximation (for short primers, NN is overkill)."""
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    if len(seq) < 14:
        return 2 * at + 4 * gc
    return 64.9 + 41 * (gc - 16.4) / (at + gc)


def _gc(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq) * 100 if seq else 0.0


def _has_hairpin(seq: str, min_stem: int = 4) -> bool:
    seq = seq.upper()
    comp = str.maketrans("ACGT", "TGCA")
    rc = seq.translate(comp)[::-1]
    for i in range(len(seq) - min_stem * 2 - 3):
        stem = seq[i:i+min_stem]
        if stem in rc:
            return True
    return False


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

    for length in range(req.primer_len_min, req.primer_len_max + 1):
        for pos in range(0, len(template) - length + 1):
            seq = template[pos:pos + length]
            gc = _gc(seq)
            tm = _tm_nearest_neighbor(seq)
            if (req.tm_min <= tm <= req.tm_max
                    and req.gc_min <= gc <= req.gc_max
                    and not _has_hairpin(seq)
                    and not seq[-3:].count("G") + seq[-3:].count("C") > 2  # 3' GC clamp ≤2
                    and "AAAA" not in seq and "TTTT" not in seq
                    and "GGGG" not in seq and "CCCC" not in seq):
                fwd_candidates.append(Primer(
                    sequence=seq, position=pos, length=length,
                    tm=round(tm, 1), gc_content=round(gc, 1), direction="forward",
                ))

        for pos in range(length, len(template) + 1):
            seq = template[pos - length:pos]
            rc_seq = _reverse_complement(seq)
            gc = _gc(rc_seq)
            tm = _tm_nearest_neighbor(rc_seq)
            if (req.tm_min <= tm <= req.tm_max
                    and req.gc_min <= gc <= req.gc_max
                    and not _has_hairpin(rc_seq)):
                rev_candidates.append(Primer(
                    sequence=rc_seq, position=pos - length, length=length,
                    tm=round(tm, 1), gc_content=round(gc, 1), direction="reverse",
                ))

    pairs: list[PrimerPair] = []
    for fwd in fwd_candidates:
        for rev in rev_candidates:
            product = rev.position + rev.length - fwd.position
            if not (req.product_min <= product <= req.product_max):
                continue
            if rev.position <= fwd.position:
                continue
            tm_diff = abs(fwd.tm - rev.tm)
            if tm_diff > 5:
                continue
            penalty = tm_diff + abs(fwd.gc_content - rev.gc_content) * 0.1
            pairs.append(PrimerPair(
                forward=fwd, reverse=rev,
                product_size=product,
                penalty=round(penalty, 3),
            ))

    pairs.sort(key=lambda p: p.penalty)
    return pairs[:req.max_pairs]


# Common restriction enzymes with recognition sequences
RESTRICTION_ENZYMES = {
    "EcoRI":  "GAATTC",
    "BamHI":  "GGATCC",
    "HindIII": "AAGCTT",
    "SalI":   "GTCGAC",
    "XbaI":   "TCTAGA",
    "SmaI":   "CCCGGG",
    "KpnI":   "GGTACC",
    "SacI":   "GAGCTC",
    "NotI":   "GCGGCCGC",
    "XhoI":   "CTCGAG",
    "NcoI":   "CCATGG",
    "NdeI":   "CATATG",
    "ClaI":   "ATCGAT",
    "SphI":   "GCATGC",
    "PstI":   "CTGCAG",
    "PvuII":  "CAGCTG",
    "AvaI":   "CYCGRG",
    "AvaII":  "GGWCC",
    "EcoRV":  "GATATC",
    "MluI":   "ACGCGT",
    "NheI":   "GCTAGC",
    "AgeI":   "ACCGGT",
    "BglII":  "AGATCT",
    "MfeI":   "CAATTG",
}

IUPAC = {"R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
         "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
         "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}


def _pattern(recog: str) -> str:
    return "".join(IUPAC.get(c, c) for c in recog.upper())


class RestrictionRequest(BaseModel):
    template: str
    enzymes: list[str] | None = None


@router.post("/restriction_sites", response_model=list[RestrictionSite])
def restriction_sites(req: RestrictionRequest) -> list[RestrictionSite]:
    template, enzymes = req.template, req.enzymes
    """Find restriction enzyme cut sites in a sequence."""
    template = template.upper().replace(" ", "").replace("\n", "")  # noqa: F841 (reassigned from req)
    target_enzymes = {k: v for k, v in RESTRICTION_ENZYMES.items()
                      if not enzymes or k in enzymes}
    results = []
    for name, recog in target_enzymes.items():
        pattern = _pattern(recog)
        positions = [m.start() for m in re.finditer(f"(?={pattern})", template)]
        if positions:
            results.append(RestrictionSite(
                enzyme=name,
                cut_pattern=recog,
                positions=positions,
                count=len(positions),
            ))
    results.sort(key=lambda r: (-r.count, r.enzyme))
    return results
