"""GenBank / EMBL format parsing and writing via Biopython."""

from io import StringIO
from typing import TextIO, BinaryIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ..models.sequence import Sequence, MoleculeType, Annotation, Strand
from .insdc_location import read_location, write_location


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
            # Every interval, in file order. A CompoundLocation reports each part;
            # a simple FeatureLocation reports itself as a single part. Reading
            # only start/end collapsed every join() to its outer bounds (#57).
            # `read_location` also keeps what the integers cannot carry: the
            # `<`/`>` of each boundary, the join-vs-order operator and a part's
            # remote accession, all of which used to be dropped here (#94).
            parts, part_details, operator = read_location(feat.location)
        except Exception:
            continue
        annotations.append(Annotation(
            feature_type=feat.type,
            start=start,
            end=end,
            parts=parts or [(start, end)],
            part_details=part_details,
            location_operator=operator,
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
    from Bio.SeqFeature import SeqFeature
    records = []
    for s in sequences:
        seq = Seq(s.seq)
        r = SeqRecord(seq, id=s.accession or s.id, name=s.name or s.id[:16], description=s.description)
        r.annotations["molecule_type"] = s.molecule_type.value
        r.annotations["topology"] = "circular" if s.is_circular else "linear"
        for ann in s.annotations:
            strand_val = 1 if ann.strand == Strand.PLUS else -1 if ann.strand == Strand.MINUS else 0

            # Rebuild the real location. More than one part means a compound
            # location; the old writer always emitted a single FeatureLocation
            # over the outer bounds, which is the export half of #57. It then
            # emitted only the coordinates, dropping every partial boundary,
            # rewriting order() as join() and localising remote parts (#94) —
            # `write_location` puts all three back.
            location = write_location(ann, strand_val)

            feat = SeqFeature(
                location,
                type=ann.feature_type,
                # A list qualifier stays a list: `[str(v)]` turned two /note
                # entries into one /note="['first', 'second']" (#57).
                qualifiers={
                    k: [str(x) for x in v] if isinstance(v, (list, tuple)) else [str(v)]
                    for k, v in ann.qualifiers.items()
                },
            )
            r.features.append(feat)
        records.append(r)
    buf = StringIO()
    SeqIO.write(records, buf, "genbank")
    return buf.getvalue()
