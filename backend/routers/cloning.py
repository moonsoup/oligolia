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
