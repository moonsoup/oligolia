"""Molecular cloning simulation — restriction-ligation assembly (issue #35).

Joins digested fragments (carrying the overhang data produced by the digest
endpoint, issue #34) into a new circular sequence, enforcing overhang
compatibility at every junction. Because the digest representation places each
overhang's bases on exactly one fragment's top strand, a valid ligation
reconstructs by ordered top-strand concatenation — this module's job is to
verify the junctions actually anneal before joining.
"""

from __future__ import annotations

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from ..models.sequence import MoleculeType, Sequence

router = APIRouter(prefix="/cloning", tags=["cloning"])

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


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


def _ends_compatible(
    up_seq: str, up_type: str, down_seq: str, down_type: str
) -> tuple[bool, str, str]:
    """Can upstream's right end ligate to downstream's left end?

    Returns (compatible, junction_kind, annealing_overhang).
    """
    if up_type == "blunt" and down_type == "blunt":
        return True, "blunt", ""
    if up_type in ("5'", "3'") and up_type == down_type:
        # The two single-stranded overhangs anneal iff they are reverse
        # complements. For this panel's palindromic-site enzymes the stored
        # top-strand bases are self-complementary, so this reduces to equality.
        if up_seq and up_seq == _reverse_complement(down_seq):
            return True, f"sticky-{up_type}", up_seq
    return False, "", ""


@router.post("/ligate", response_model=LigationResult)
def ligate(req: LigationRequest) -> LigationResult:
    """Ligate fragments in the given order, checking every junction."""
    frags = req.fragments

    # Adjacent junctions, plus the closing junction (last -> first) if circular.
    pairs = [(i, i + 1) for i in range(len(frags) - 1)]
    if req.circular:
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

    product_seq = "".join(f.sequence for f in frags)
    product = Sequence(
        id=req.product_name,
        name=req.product_name,
        seq=product_seq,
        molecule_type=MoleculeType.DNA,
        is_circular=req.circular,
    )
    return LigationResult(product=product, junctions=junctions)


# ── Gibson assembly (issue #36) ───────────────────────────────────────────────

# Exhaustive ordering search is factorial; refuse absurd fragment counts.
_MAX_GIBSON_FRAGMENTS = 12


class GibsonRequest(BaseModel):
    fragments: list[str] = Field(min_length=2)
    min_overlap: int = Field(default=15, ge=5)
    product_name: str = "gibson_product"


class GibsonJunction(BaseModel):
    upstream_index: int    # index into the input fragment list
    downstream_index: int
    overlap_length: int
    overlap: str


class GibsonResult(BaseModel):
    product: Sequence
    order: list[int]                 # input indices in assembled circular order
    junctions: list[GibsonJunction]


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
    frags = [f.upper().replace(" ", "").replace("\n", "") for f in req.fragments]
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

    seq = Sequence(
        id=req.product_name, name=req.product_name, seq=product,
        molecule_type=MoleculeType.DNA, is_circular=True,
    )
    return GibsonResult(product=seq, order=order, junctions=junctions)
