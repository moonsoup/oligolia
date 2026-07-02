"""Molecular cloning simulation — restriction-ligation assembly (issue #35).

Joins digested fragments (carrying the overhang data produced by the digest
endpoint, issue #34) into a new circular sequence, enforcing overhang
compatibility at every junction. Because the digest representation places each
overhang's bases on exactly one fragment's top strand, a valid ligation
reconstructs by ordered top-strand concatenation — this module's job is to
verify the junctions actually anneal before joining.
"""

from __future__ import annotations

from Bio import Restriction
from Bio.Seq import Seq
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field, field_validator

from ..models.sequence import Annotation, MoleculeType, Sequence
from .primers import _overhang_at

router = APIRouter(prefix="/cloning", tags=["cloning"])

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _reindex_annotations(
    annotations: list[Annotation],
    local_lo: int,
    local_hi: int,
    product_offset: int,
    label: str,
    warnings: list[str],
) -> list[Annotation]:
    """Map a fragment's annotations onto product coordinates.

    Only the fragment's retained window ``[local_lo, local_hi)`` (its bases
    that actually survive into the product — trimmed overlaps are excluded)
    is placed, starting at ``product_offset``. Junction-split policy: an
    annotation extending past the retained window is **truncated** to the
    retained portion and flagged ``qualifiers["truncated_at_junction"]="true"``;
    one lying entirely outside the window is **dropped**. Both cases emit a
    warning so nothing is silently altered.
    """
    out: list[Annotation] = []
    for ann in annotations:
        clip_s, clip_e = max(ann.start, local_lo), min(ann.end, local_hi)
        if clip_s >= clip_e:
            warnings.append(
                f"{label}: dropped {ann.feature_type} [{ann.start}:{ann.end}] "
                "(outside retained region / in a trimmed overlap)"
            )
            continue
        new_s = product_offset + (clip_s - local_lo)
        new_e = product_offset + (clip_e - local_lo)
        qualifiers = dict(ann.qualifiers)
        if clip_s != ann.start or clip_e != ann.end:
            qualifiers["truncated_at_junction"] = "true"
            warnings.append(
                f"{label}: truncated {ann.feature_type} [{ann.start}:{ann.end}] "
                f"to [{new_s}:{new_e}] at an assembly junction"
            )
        out.append(Annotation(
            feature_type=ann.feature_type, start=new_s, end=new_e,
            strand=ann.strand, qualifiers=qualifiers,
        ))
    return out


class LigationFragment(BaseModel):
    """A digested fragment ready for ligation.

    Field names/shape match ``DigestFragment`` so digest output can be fed in
    directly. ``*_overhang_type`` is one of "5'", "3'", "blunt", or "none"
    (a free linear terminus, which cannot form a defined junction).
    """

    sequence: str
    name: str = ""
    left_overhang: str = ""
    left_overhang_type: str = "none"
    right_overhang: str = ""
    right_overhang_type: str = "none"
    # Source annotations in this fragment's local coordinates; reindexed onto
    # the product (issue #38).
    annotations: list[Annotation] = Field(default_factory=list)


class LigationRequest(BaseModel):
    fragments: list[LigationFragment] = Field(min_length=2)
    circular: bool = True  # ligate into a circle (also checks last->first)
    product_name: str = "ligation_product"


class Junction(BaseModel):
    upstream: str
    downstream: str
    kind: str          # "sticky-5'", "sticky-3'", or "blunt"
    overhang: str      # annealing overhang bases (empty for blunt)


class LigationResult(BaseModel):
    product: Sequence
    junctions: list[Junction]
    warnings: list[str] = Field(default_factory=list)


def _ends_compatible(
    up_seq: str, up_type: str, down_seq: str, down_type: str
) -> tuple[bool, str, str]:
    """Can upstream's right end ligate to downstream's left end?

    Returns (compatible, junction_kind, annealing_overhang).
    """
    if up_type == "blunt" and down_type == "blunt":
        return True, "blunt", ""
    if up_type in ("5'", "3'") and up_type == down_type:
        # Overhangs are stored as top-strand stagger bases, so a compatible
        # junction is one where both ends carry the same 4 bp — i.e. the same
        # sequence reads across the junction on the top strand. This also
        # correctly matches isoschizomer-compatible ends (e.g. BamHI/BglII,
        # both GATC) and Golden Gate's designed non-palindromic overhangs.
        if up_seq and up_seq == down_seq:
            return True, f"sticky-{up_type}", up_seq
    return False, "", ""


def _ligate_fragments(
    frags: list[LigationFragment], circular: bool, product_name: str
) -> LigationResult:
    """Core assembly: check junctions in order, concatenate, reindex annotations.

    Shared by restriction-ligation (/ligate) and Golden Gate (/goldengate).
    Raises 400 on the first incompatible junction.
    """
    # Adjacent junctions, plus the closing junction (last -> first) if circular.
    pairs = [(i, i + 1) for i in range(len(frags) - 1)]
    if circular:
        pairs.append((len(frags) - 1, 0))

    junctions: list[Junction] = []
    for up_i, down_i in pairs:
        up, down = frags[up_i], frags[down_i]
        ok, kind, overhang = _ends_compatible(
            up.right_overhang, up.right_overhang_type,
            down.left_overhang, down.left_overhang_type,
        )
        up_name = up.name or f"fragment {up_i + 1}"
        down_name = down.name or f"fragment {down_i + 1}"
        if not ok:
            raise HTTPException(
                400,
                f"Incompatible ends between {up_name} (3' end: "
                f"{up.right_overhang_type} '{up.right_overhang}') and {down_name} "
                f"(5' end: {down.left_overhang_type} '{down.left_overhang}'). "
                "Fragments cannot be ligated as ordered.",
            )
        junctions.append(Junction(
            upstream=up_name, downstream=down_name, kind=kind, overhang=overhang,
        ))

    # Reindex source annotations onto product coordinates. Ligation
    # concatenates full top strands, so each fragment's placed window is its
    # whole length at the running offset — no overlap trimming, no splits.
    warnings: list[str] = []
    annotations: list[Annotation] = []
    offset = 0
    for i, f in enumerate(frags):
        label = f.name or f"fragment {i + 1}"
        annotations += _reindex_annotations(
            f.annotations, 0, len(f.sequence), offset, label, warnings)
        offset += len(f.sequence)

    product = Sequence(
        id=product_name, name=product_name,
        seq="".join(f.sequence for f in frags),
        molecule_type=MoleculeType.DNA, is_circular=circular,
        annotations=annotations,
    )
    return LigationResult(product=product, junctions=junctions, warnings=warnings)


@router.post("/ligate", response_model=LigationResult)
def ligate(req: LigationRequest) -> LigationResult:
    """Ligate fragments in the given order, checking every junction."""
    return _ligate_fragments(req.fragments, req.circular, req.product_name)


# ── Gibson assembly (issue #36) ───────────────────────────────────────────────

# Exhaustive ordering search is factorial; refuse absurd fragment counts.
_MAX_GIBSON_FRAGMENTS = 12


class GibsonFragment(BaseModel):
    sequence: str
    name: str = ""
    annotations: list[Annotation] = Field(default_factory=list)


class GibsonRequest(BaseModel):
    fragments: list[GibsonFragment] = Field(min_length=2)
    min_overlap: int = Field(default=15, ge=5)
    product_name: str = "gibson_product"

    @field_validator("fragments", mode="before")
    @classmethod
    def _coerce_bare_strings(cls, v: object) -> object:
        # Accept a plain list of sequence strings for convenience/back-compat.
        if isinstance(v, list):
            return [{"sequence": x} if isinstance(x, str) else x for x in v]
        return v


class GibsonJunction(BaseModel):
    upstream_index: int    # index into the input fragment list
    downstream_index: int
    overlap_length: int
    overlap: str


class GibsonResult(BaseModel):
    product: Sequence
    order: list[int]                 # input indices in assembled circular order
    junctions: list[GibsonJunction]
    warnings: list[str] = Field(default_factory=list)


def _overlap(a: str, b: str, min_overlap: int) -> int:
    """Longest L in [min_overlap, min(len)] with suffix(a, L) == prefix(b, L)."""
    for length in range(min(len(a), len(b)), min_overlap - 1, -1):
        if a[-length:] == b[:length]:
            return length
    return 0


@router.post("/gibson", response_model=GibsonResult)
def gibson(req: GibsonRequest) -> GibsonResult:
    """Assemble linear fragments with overlapping ends into a circular product.

    v1 matches ends exactly (suffix→prefix) in the given orientation and
    searches all circular orderings; ambiguity or no valid ordering is
    reported as a clear error rather than guessed. Reverse-complement
    orientation and mismatch-tolerant overlaps are noted follow-ups.
    """
    gfrags = req.fragments
    frags = [g.sequence.upper().replace(" ", "").replace("\n", "") for g in gfrags]
    n = len(frags)
    if n > _MAX_GIBSON_FRAGMENTS:
        raise HTTPException(400, f"Too many fragments for v1 assembly (max {_MAX_GIBSON_FRAGMENTS}).")

    # ov[i][j] = overlap length of fragment i's suffix into fragment j's prefix.
    ov = [[_overlap(frags[i], frags[j], req.min_overlap) if i != j else 0
           for j in range(n)] for i in range(n)]

    # Find every circular ordering (anchored at fragment 0 to factor out
    # rotation) where consecutive fragments overlap and the circle closes.
    orders: list[list[int]] = []

    def _dfs(path: list[int], used: set[int]) -> None:
        if len(path) == n:
            if ov[path[-1]][path[0]] > 0:
                orders.append(list(path))
            return
        for nxt in range(n):
            if nxt not in used and ov[path[-1]][nxt] > 0:
                path.append(nxt)
                used.add(nxt)
                _dfs(path, used)
                used.discard(nxt)
                path.pop()

    _dfs([0], {0})

    if not orders:
        raise HTTPException(
            400,
            f"No valid circular assembly: fragments do not form a closed loop "
            f"with overlaps ≥ {req.min_overlap} bp. Check overlap design/orientation.",
        )
    if len(orders) > 1:
        raise HTTPException(
            400,
            f"Ambiguous assembly: {len(orders)} distinct circular orderings satisfy "
            "the overlaps. Refusing to guess — make overlap junctions unique.",
        )

    order = orders[0]
    junctions: list[GibsonJunction] = []
    product = frags[order[0]]
    for k in range(1, n):
        length = ov[order[k - 1]][order[k]]
        junctions.append(GibsonJunction(
            upstream_index=order[k - 1], downstream_index=order[k],
            overlap_length=length, overlap=frags[order[k]][:length],
        ))
    # Rebuild trimming leading overlaps; then trim the closing overlap once so
    # the shared region isn't duplicated at both ends of the linear form.
    for k in range(1, n):
        product += frags[order[k]][ov[order[k - 1]][order[k]]:]
    closing = ov[order[-1]][order[0]]
    junctions.append(GibsonJunction(
        upstream_index=order[-1], downstream_index=order[0],
        overlap_length=closing, overlap=frags[order[0]][:closing],
    ))
    if closing:
        product = product[:-closing]

    # Reindex source annotations onto product coordinates. Each fragment's
    # retained window excludes its leading overlap (already contributed by the
    # previous fragment) and, for the last fragment, the trailing closing
    # overlap (duplicated by the first fragment's prefix). See #38 policy.
    warnings: list[str] = []
    annotations: list[Annotation] = []
    running = 0
    for k in range(n):
        idx = order[k]
        local_lo = 0 if k == 0 else ov[order[k - 1]][order[k]]
        local_hi = len(frags[idx]) - (closing if k == n - 1 else 0)
        label = gfrags[idx].name or f"fragment {idx + 1}"
        annotations += _reindex_annotations(
            gfrags[idx].annotations, local_lo, local_hi, running, label, warnings)
        running += local_hi - local_lo

    seq = Sequence(
        id=req.product_name, name=req.product_name, seq=product,
        molecule_type=MoleculeType.DNA, is_circular=True,
        annotations=annotations,
    )
    return GibsonResult(product=seq, order=order, junctions=junctions, warnings=warnings)


# ── Golden Gate assembly (issue #37) ──────────────────────────────────────────

# Type IIS enzymes cut outside their recognition site, leaving a user-defined
# 4 bp overhang — the basis of scarless, ordered Golden Gate assembly. Cut
# geometry comes from Bio.Restriction (verified to model the remote cut).
_TYPE_IIS = {
    name: getattr(Restriction, name)
    for name in ("BsaI", "BsmBI", "Esp3I", "BbsI")
    if hasattr(Restriction, name)
}


class GoldenGateRequest(BaseModel):
    # Parts in junction order; each carries two inward-facing Type IIS sites.
    parts: list[str] = Field(min_length=2)
    enzyme: str = "BsaI"
    product_name: str = "goldengate_product"


def _digest_linear_fragments(seq: str, enz) -> list[dict]:
    """Linear digest of ``seq`` with a Bio.Restriction enzyme, with overhangs.

    Mirrors the digest geometry in routers/primers so the released Golden Gate
    insert carries the same overhang representation the ligation core expects.
    """
    cuts: dict[int, tuple[str, str]] = {}
    for pos in enz.search(Seq(seq), linear=True):
        c = pos - 1
        cuts[c] = _overhang_at(seq, c, enz.ovhg, False)
    boundaries = [0, *sorted(cuts), len(seq)]
    frags = []
    for i in range(len(boundaries) - 1):
        s, e = boundaries[i], boundaries[i + 1]
        if e <= s:
            continue
        lo, lt = cuts.get(s, ("", "none"))
        ro, rt = cuts.get(e, ("", "none"))
        frags.append({
            "sequence": seq[s:e],
            "left_overhang": lo, "left_overhang_type": lt,
            "right_overhang": ro, "right_overhang_type": rt,
        })
    return frags


@router.post("/goldengate", response_model=LigationResult)
def golden_gate(req: GoldenGateRequest) -> LigationResult:
    """Simulate Golden Gate: Type IIS-release each part's insert, then ligate.

    Each part is digested with the chosen Type IIS enzyme; the single fragment
    cut on both ends (recognition sites removed) is the released insert. Inserts
    are ligated in the given part order into a circular product via the shared
    overhang-checked assembly. Junction overhangs are validated to be unique and
    non-palindromic (standard Golden Gate constraints) and flagged otherwise.
    """
    enz = _TYPE_IIS.get(req.enzyme)
    if enz is None:
        raise HTTPException(
            400, f"'{req.enzyme}' is not a supported Type IIS enzyme "
            f"(choose from {sorted(_TYPE_IIS)}).")

    inserts: list[LigationFragment] = []
    for i, raw in enumerate(req.parts):
        part = raw.upper().replace(" ", "").replace("\n", "")
        frags = _digest_linear_fragments(part, enz)
        internal = [f for f in frags
                    if f["left_overhang_type"] != "none" and f["right_overhang_type"] != "none"]
        if len(internal) != 1:
            raise HTTPException(
                400,
                f"Part {i + 1} released {len(internal)} inserts (expected exactly 1). "
                f"A Golden Gate part needs two inward-facing {req.enzyme} sites.",
            )
        ins = internal[0]
        inserts.append(LigationFragment(name=f"insert {i + 1}", **ins))

    # Golden Gate design constraints: junction overhangs unique + non-palindromic.
    warnings: list[str] = []
    junction_overhangs = [ins.right_overhang for ins in inserts]
    dupes = sorted({o for o in junction_overhangs if junction_overhangs.count(o) > 1})
    if dupes:
        warnings.append(
            f"Colliding junction overhangs {dupes}: Golden Gate junctions must be "
            "unique or fragments can assemble in the wrong order.")
    for o in junction_overhangs:
        if o and o == _reverse_complement(o):
            warnings.append(
                f"Palindromic junction overhang '{o}': can ligate in either "
                "orientation, making assembly ambiguous.")

    result = _ligate_fragments(inserts, circular=True, product_name=req.product_name)
    result.warnings = warnings + result.warnings
    return result
