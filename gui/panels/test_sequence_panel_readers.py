"""#67.2: `.embl` was routed to read_genbank, so EMBL files loaded as nothing.

`read_embl` existed and worked the whole time. Biopython's GenBank parser simply
yields no records for an EMBL file, and `for seq in seqs:` over an empty list is
silent — so the user saw a file open and nothing appear.

The reader table is module-level so the routing is testable without driving a
file dialog.
"""

from __future__ import annotations

import os
from io import StringIO

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from backend.formats import read_embl, read_fasta, read_fastq, read_genbank, read_snapgene  # noqa: E402
from gui.panels.sequence_panel import BINARY_EXTENSIONS, READERS_BY_EXTENSION  # noqa: E402

EMBL = """ID   DEMO; SV 1; linear; genomic DNA; STD; SYN; 60 BP.
AC   DEMO;
DE   a small embl record
FH   Key             Location/Qualifiers
FH
FT   source          1..60
FT   CDS             5..40
FT                   /gene="demo"
SQ   Sequence 60 BP; 15 A; 15 C; 15 G; 15 T; 0 other;
     atggcctgtg ggcatttggc caatttaggc catggacgtg gcatcacgtg gcatcacgtt        60
//
"""


def test_embl_is_routed_to_the_embl_reader() -> None:
    assert READERS_BY_EXTENSION["embl"] is read_embl


def test_every_other_extension_keeps_its_reader() -> None:
    """Guard against fixing one route by breaking another."""
    assert READERS_BY_EXTENSION["fasta"] is read_fasta
    assert READERS_BY_EXTENSION["fa"] is read_fasta
    assert READERS_BY_EXTENSION["fna"] is read_fasta
    assert READERS_BY_EXTENSION["faa"] is read_fasta
    assert READERS_BY_EXTENSION["fastq"] is read_fastq
    assert READERS_BY_EXTENSION["fq"] is read_fastq
    assert READERS_BY_EXTENSION["gb"] is read_genbank
    assert READERS_BY_EXTENSION["gbk"] is read_genbank
    assert READERS_BY_EXTENSION["dna"] is read_snapgene


def test_only_snapgene_is_read_in_binary() -> None:
    assert BINARY_EXTENSIONS == {"dna"}


def test_the_embl_reader_actually_reads_this_record() -> None:
    """Routing is only half of it — the reader has to work."""
    seqs = read_embl(StringIO(EMBL))
    assert len(seqs) == 1, seqs
    seq = seqs[0]
    assert len(seq.seq) == 60, seq.seq
    assert any(a.feature_type == "CDS" for a in seq.annotations), seq.annotations


def test_the_genbank_reader_really_does_return_nothing_for_embl() -> None:
    """The reason this was silent: no exception, just an empty list."""
    assert read_genbank(StringIO(EMBL)) == []


def test_the_embl_cds_keeps_its_location_parts() -> None:
    """#57's parts field should hold through the EMBL path too."""
    cds = next(a for a in read_embl(StringIO(EMBL))[0].annotations if a.feature_type == "CDS")
    assert cds.parts == [(4, 40)], cds.parts
