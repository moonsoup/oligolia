"""Tests for SnapGene .dna import (issue #41).

Biopython has no .dna writer, so we synthesize a minimal-but-valid SnapGene
binary from its documented packet structure (cookie / DNA / Features) and
read it back through our reader.
"""

import struct

from ..formats import read_snapgene
from ..models.sequence import MoleculeType

_FEATURES_XML = (
    '<Features>'
    '<Feature type="CDS" name="testGene" directionality="1">'
    '<Segment range="1-9" type="standard"/>'
    '<Q name="label"><V text="testGene"/></Q>'
    '</Feature>'
    '</Features>'
)


def _packet(ptype: int, data: bytes) -> bytes:
    return struct.pack(">B", ptype) + struct.pack(">I", len(data)) + data


def _make_dna(seq: str, circular: bool, features_xml: str | None = None) -> bytes:
    cookie = _packet(0x09, struct.pack(">8sHHH", b"SnapGene", 1, 1, 1))
    dna = _packet(0x00, struct.pack(">B", 1 if circular else 0) + seq.encode("ascii"))
    out = cookie + dna
    if features_xml:
        out += _packet(0x0A, features_xml.encode("utf-8"))
    return out


def test_reads_circular_dna_with_features() -> None:
    seq = "ATGGCCTAAGGGCCCTTTAAA"  # 21 bp
    data = _make_dna(seq, circular=True, features_xml=_FEATURES_XML)
    seqs = read_snapgene(data)
    assert len(seqs) == 1
    s = seqs[0]
    assert s.seq == seq
    assert s.molecule_type == MoleculeType.DNA
    assert s.is_circular is True  # topology carried through (issue #40 reuse)
    cds = [a for a in s.annotations if a.feature_type == "CDS"]
    assert len(cds) == 1
    assert cds[0].start == 0 and cds[0].end == 9  # SnapGene 1-based -> 0-based
    assert cds[0].qualifiers.get("label") == "testGene"


def test_reads_linear_dna() -> None:
    data = _make_dna("ACGTACGTACGT", circular=False)
    s = read_snapgene(data)[0]
    assert s.is_circular is False
    assert s.seq == "ACGTACGTACGT"


def test_reads_from_file_path(tmp_path) -> None:
    p = tmp_path / "plasmid.dna"
    p.write_bytes(_make_dna("GGGGCCCCAAAATTTT", circular=True))
    s = read_snapgene(str(p))[0]
    assert s.is_circular is True
    assert s.seq == "GGGGCCCCAAAATTTT"


def test_reads_from_binary_handle(tmp_path) -> None:
    p = tmp_path / "plasmid.dna"
    p.write_bytes(_make_dna("ATATATATGCGC", circular=False))
    with open(p, "rb") as fh:
        s = read_snapgene(fh)[0]
    assert s.seq == "ATATATATGCGC"
