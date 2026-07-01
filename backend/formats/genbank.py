"""GenBank / EMBL format parsing and writing via Biopython."""

from io import StringIO
from typing import TextIO, BinaryIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ..models.sequence import Sequence, MoleculeType, Annotation, Strand


def _strand(val: int | None) -> Strand:
    if val == 1:
        return Strand.PLUS
    if val == -1:
        return Strand.MINUS
    return Strand.BOTH


def _to_sequence(r: SeqRecord) -> Sequence:
    mol_raw = str(r.annotations.get("molecule_type", "")).upper()
    if "DNA" in mol_raw:
        mol = MoleculeType.DNA
    elif "RNA" in mol_raw:
        mol = MoleculeType.RNA
    elif "PROTEIN" in mol_raw or "AA" in mol_raw:
        mol = MoleculeType.PROTEIN
    else:
        mol = MoleculeType.UNKNOWN

    annotations = []
    for feat in r.features:
        try:
            start = int(feat.location.start)
            end = int(feat.location.end)
        except Exception:
            continue
        annotations.append(Annotation(
            feature_type=feat.type,
            start=start,
            end=end,
            strand=_strand(feat.location.strand),
            qualifiers={k: (v[0] if isinstance(v, list) and len(v) == 1 else v)
                        for k, v in feat.qualifiers.items()},
        ))

    return Sequence(
        id=r.id,
        name=r.name,
        description=r.description,
        seq=str(r.seq),
        molecule_type=mol,
        annotations=annotations,
        accession=r.id,
        source_db="genbank",
        is_circular=r.annotations.get("topology") == "circular",
    )


def read_genbank(source: str | TextIO | BinaryIO) -> list[Sequence]:
    if isinstance(source, str):
        records = list(SeqIO.parse(StringIO(source), "genbank"))
    else:
        records = list(SeqIO.parse(source, "genbank"))
    return [_to_sequence(r) for r in records]


def read_embl(source: str | TextIO | BinaryIO) -> list[Sequence]:
    if isinstance(source, str):
        records = list(SeqIO.parse(StringIO(source), "embl"))
    else:
        records = list(SeqIO.parse(source, "embl"))
    return [_to_sequence(r) for r in records]


def write_genbank(sequences: list[Sequence]) -> str:
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    records = []
    for s in sequences:
        seq = Seq(s.seq)
        r = SeqRecord(seq, id=s.accession or s.id, name=s.name or s.id[:16], description=s.description)
        r.annotations["molecule_type"] = s.molecule_type.value
        r.annotations["topology"] = "circular" if s.is_circular else "linear"
        for ann in s.annotations:
            strand_val = 1 if ann.strand == Strand.PLUS else -1 if ann.strand == Strand.MINUS else 0
            feat = SeqFeature(
                FeatureLocation(ann.start, ann.end, strand=strand_val),
                type=ann.feature_type,
                qualifiers={k: [str(v)] for k, v in ann.qualifiers.items()},
            )
            r.features.append(feat)
        records.append(r)
    buf = StringIO()
    SeqIO.write(records, buf, "genbank")
    return buf.getvalue()
