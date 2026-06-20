"""FASTA/FASTQ format read/write via Biopython."""

from io import StringIO
from typing import BinaryIO, TextIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ..models.sequence import Sequence, MoleculeType, Annotation


def _molecule_type(seq: SeqRecord) -> MoleculeType:
    mol = str(seq.annotations.get("molecule_type", "")).upper()
    if "DNA" in mol:
        return MoleculeType.DNA
    if "RNA" in mol:
        return MoleculeType.RNA
    if "PROTEIN" in mol or "AA" in mol:
        return MoleculeType.PROTEIN
    # Infer from alphabet characters
    upper = str(seq.seq).upper()
    chars = set(upper) - {"-", "N", "*", "X"}
    rna_only = chars - set("ACGT")
    if "U" in chars and "T" not in chars:
        return MoleculeType.RNA
    if rna_only - {"U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}:
        return MoleculeType.PROTEIN
    return MoleculeType.DNA


def read_fasta(source: str | TextIO | BinaryIO) -> list[Sequence]:
    if isinstance(source, str):
        records = list(SeqIO.parse(StringIO(source), "fasta"))
    else:
        records = list(SeqIO.parse(source, "fasta"))
    return [
        Sequence(
            id=r.id,
            name=r.name,
            description=r.description,
            seq=str(r.seq),
            molecule_type=_molecule_type(r),
        )
        for r in records
    ]


def read_fastq(source: str | TextIO | BinaryIO) -> list[Sequence]:
    if isinstance(source, str):
        records = list(SeqIO.parse(StringIO(source), "fastq"))
    else:
        records = list(SeqIO.parse(source, "fastq"))
    return [
        Sequence(
            id=r.id,
            name=r.name,
            description=r.description,
            seq=str(r.seq),
            molecule_type=MoleculeType.DNA,
            annotations=[
                Annotation(
                    feature_type="quality",
                    start=0,
                    end=len(r.seq),
                    qualifiers={"phred_quality": r.letter_annotations.get("phred_quality", [])},
                )
            ],
        )
        for r in records
    ]


def write_fasta(sequences: list[Sequence]) -> str:
    buf = StringIO()
    for s in sequences:
        buf.write(f">{s.id} {s.description}\n")
        seq = s.seq
        for i in range(0, len(seq), 60):
            buf.write(seq[i:i+60] + "\n")
    return buf.getvalue()


def write_fastq(sequences: list[Sequence], default_quality: int = 40) -> str:
    """Write FASTQ; quality from annotation if present, else uniform default_quality."""
    buf = StringIO()
    for s in sequences:
        quality_scores: list[int] = []
        for ann in s.annotations:
            if ann.feature_type == "quality":
                quality_scores = ann.qualifiers.get("phred_quality", [])
                break
        if not quality_scores:
            quality_scores = [default_quality] * len(s.seq)
        buf.write(f"@{s.id} {s.description}\n")
        buf.write(s.seq + "\n")
        buf.write("+\n")
        buf.write("".join(chr(q + 33) for q in quality_scores) + "\n")
    return buf.getvalue()
