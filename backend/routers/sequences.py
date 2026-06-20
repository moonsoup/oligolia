"""Sequence editing and manipulation endpoints."""

from fastapi import APIRouter, HTTPException
from Bio.Seq import Seq
from ..models.sequence import Sequence, SequenceEditRequest, SequenceEditResult, MoleculeType

router = APIRouter(prefix="/sequences", tags=["sequences"])

# In-memory session store (keyed by id)
_store: dict[str, Sequence] = {}


@router.get("/", response_model=list[Sequence])
def list_sequences() -> list[Sequence]:
    return list(_store.values())


@router.get("/{seq_id}", response_model=Sequence)
def get_sequence(seq_id: str) -> Sequence:
    if seq_id not in _store:
        raise HTTPException(404, f"Sequence {seq_id!r} not found")
    return _store[seq_id]


@router.post("/", response_model=Sequence, status_code=201)
def add_sequence(seq: Sequence) -> Sequence:
    _store[seq.id] = seq
    return seq


@router.delete("/{seq_id}", status_code=204)
def delete_sequence(seq_id: str) -> None:
    if seq_id not in _store:
        raise HTTPException(404, f"Sequence {seq_id!r} not found")
    del _store[seq_id]


@router.post("/{seq_id}/edit", response_model=SequenceEditResult)
def edit_sequence(seq_id: str, req: SequenceEditRequest) -> SequenceEditResult:
    if seq_id not in _store:
        raise HTTPException(404, f"Sequence {seq_id!r} not found")
    s = _store[seq_id]
    seq_str = s.seq
    op = req.operation.lower()

    if op == "insert":
        if req.position is None or req.insert_seq is None:
            raise HTTPException(400, "insert requires position and insert_seq")
        pos = req.position
        if not (0 <= pos <= len(seq_str)):
            raise HTTPException(400, f"position {pos} out of range [0, {len(seq_str)}]")
        new_seq = seq_str[:pos] + req.insert_seq + seq_str[pos:]
        result = SequenceEditResult(
            original_id=seq_id, operation=op, result_seq=new_seq,
            diff_start=pos, diff_end=pos + len(req.insert_seq),
            message=f"Inserted {len(req.insert_seq)} bases at position {pos}",
        )

    elif op == "delete":
        if req.position is None or req.end_position is None:
            raise HTTPException(400, "delete requires position and end_position")
        start, end = req.position, req.end_position
        if not (0 <= start < end <= len(seq_str)):
            raise HTTPException(400, f"invalid range [{start}, {end})")
        new_seq = seq_str[:start] + seq_str[end:]
        result = SequenceEditResult(
            original_id=seq_id, operation=op, result_seq=new_seq,
            diff_start=start, diff_end=start,
            message=f"Deleted {end - start} bases from position {start} to {end}",
        )

    elif op == "replace":
        if req.position is None or req.end_position is None or req.replacement is None:
            raise HTTPException(400, "replace requires position, end_position, and replacement")
        start, end = req.position, req.end_position
        if not (0 <= start <= end <= len(seq_str)):
            raise HTTPException(400, f"invalid range [{start}, {end})")
        new_seq = seq_str[:start] + req.replacement + seq_str[end:]
        result = SequenceEditResult(
            original_id=seq_id, operation=op, result_seq=new_seq,
            diff_start=start, diff_end=start + len(req.replacement),
            message=f"Replaced [{start}:{end}] with {len(req.replacement)}-base sequence",
        )

    elif op == "reverse_complement":
        bio_seq = Seq(seq_str)
        new_seq = str(bio_seq.reverse_complement())
        result = SequenceEditResult(
            original_id=seq_id, operation=op, result_seq=new_seq,
            message="Reverse complement computed",
        )

    elif op == "complement":
        bio_seq = Seq(seq_str)
        new_seq = str(bio_seq.complement())
        result = SequenceEditResult(
            original_id=seq_id, operation=op, result_seq=new_seq,
            message="Complement computed",
        )

    elif op == "translate":
        if s.molecule_type not in (MoleculeType.DNA, MoleculeType.RNA):
            raise HTTPException(400, "translate only applies to DNA/RNA sequences")
        bio_seq = Seq(seq_str)
        new_seq = str(bio_seq.translate(to_stop=True))
        result = SequenceEditResult(
            original_id=seq_id, operation=op, result_seq=new_seq,
            message=f"Translated {len(seq_str)} nt → {len(new_seq)} aa",
        )

    elif op == "transcribe":
        if s.molecule_type != MoleculeType.DNA:
            raise HTTPException(400, "transcribe only applies to DNA sequences")
        bio_seq = Seq(seq_str)
        new_seq = str(bio_seq.transcribe())
        result = SequenceEditResult(
            original_id=seq_id, operation=op, result_seq=new_seq,
            message="Transcribed DNA → RNA",
        )

    elif op == "back_transcribe":
        if s.molecule_type != MoleculeType.RNA:
            raise HTTPException(400, "back_transcribe only applies to RNA sequences")
        bio_seq = Seq(seq_str)
        new_seq = str(bio_seq.back_transcribe())
        result = SequenceEditResult(
            original_id=seq_id, operation=op, result_seq=new_seq,
            message="Back-transcribed RNA → DNA",
        )

    else:
        raise HTTPException(400, f"Unknown operation: {op!r}. Valid: insert, delete, replace, "
                                 "reverse_complement, complement, translate, transcribe, back_transcribe")

    # Persist the edited sequence under a new ID
    new_id = f"{seq_id}_{op}"
    _store[new_id] = Sequence(
        id=new_id,
        name=f"{s.name} ({op})",
        description=s.description,
        seq=result.result_seq,
        molecule_type=MoleculeType.PROTEIN if op == "translate" else
                      MoleculeType.RNA if op == "transcribe" else
                      MoleculeType.DNA if op == "back_transcribe" else s.molecule_type,
    )
    return result


@router.get("/{seq_id}/gc_content")
def gc_content(seq_id: str) -> dict:
    if seq_id not in _store:
        raise HTTPException(404, f"Sequence {seq_id!r} not found")
    seq = _store[seq_id].seq.upper()
    gc = (seq.count("G") + seq.count("C")) / len(seq) * 100 if seq else 0.0
    return {"seq_id": seq_id, "length": len(seq), "gc_content": round(gc, 2)}


@router.get("/{seq_id}/codon_usage")
def codon_usage(seq_id: str) -> dict:
    if seq_id not in _store:
        raise HTTPException(404, f"Sequence {seq_id!r} not found")
    seq = _store[seq_id].seq.upper()
    counts: dict[str, int] = {}
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3:
            counts[codon] = counts.get(codon, 0) + 1
    total = sum(counts.values())
    usage = {codon: {"count": cnt, "fraction": round(cnt / total, 4)}
             for codon, cnt in sorted(counts.items())}
    return {"seq_id": seq_id, "total_codons": total, "codon_usage": usage}


@router.get("/{seq_id}/find_motif")
def find_motif(seq_id: str, motif: str) -> dict:
    """Find all occurrences of a motif (supports IUPAC ambiguity codes)."""
    import re
    if seq_id not in _store:
        raise HTTPException(404, f"Sequence {seq_id!r} not found")
    iupac = {
        "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
        "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
        "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
    }
    pattern = "".join(iupac.get(c.upper(), re.escape(c)) for c in motif.upper())
    seq = _store[seq_id].seq.upper()
    positions = [m.start() for m in re.finditer(f"(?={pattern})", seq)]
    return {"seq_id": seq_id, "motif": motif, "pattern": pattern,
            "count": len(positions), "positions": positions}
