"""Tests for file format parsers using real biological file content."""

from ..formats import read_fasta, read_fastq, read_genbank, write_fasta
from ..formats import parse_vcf, write_vcf, parse_gff3, write_gff3
from ..models.sequence import MoleculeType


def test_read_fasta_brca2(brca2_fasta: str) -> None:
    seqs = read_fasta(brca2_fasta)
    assert len(seqs) == 1
    s = seqs[0]
    assert "NM_000059" in s.id
    assert s.seq.startswith("ATGCCTATTGG")
    assert len(s.seq) > 50


def test_read_fasta_protein(tp53_fasta: str) -> None:
    seqs = read_fasta(tp53_fasta)
    assert len(seqs) == 1
    assert seqs[0].seq.startswith("MEEPQSDP")


def test_read_fasta_roundtrip(brca2_fasta: str) -> None:
    seqs = read_fasta(brca2_fasta)
    fasta_out = write_fasta(seqs)
    seqs2 = read_fasta(fasta_out)
    assert seqs[0].seq == seqs2[0].seq


def test_read_fasta_multiple() -> None:
    multi = ">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTTT\n"
    seqs = read_fasta(multi)
    assert len(seqs) == 3
    assert seqs[0].id == "seq1"
    assert seqs[2].seq == "TTTT"


def test_read_fastq() -> None:
    fastq = "@read1\nATCGATCG\n+\nIIIIIIII\n@read2\nGCGCGCGC\n+\nAAAAAAAA\n"
    seqs = read_fastq(fastq)
    assert len(seqs) == 2
    assert seqs[0].seq == "ATCGATCG"
    assert seqs[0].molecule_type == MoleculeType.DNA
    # Quality scores
    quality = seqs[0].annotations[0].qualifiers.get("phred_quality", [])
    assert quality[0] == 40  # 'I' = ASCII 73 - 33 = 40


def test_read_genbank_real() -> None:
    # Minimal GenBank flat file format with real annotation structure
    gb_content = """LOCUS       NM_000059               140 bp    mRNA    linear   PRI 01-JAN-2024
DEFINITION  Homo sapiens BRCA2 mRNA (partial).
ACCESSION   NM_000059
VERSION     NM_000059.4
FEATURES             Location/Qualifiers
     gene            1..140
                     /gene="BRCA2"
                     /gene_synonym="FANCD1"
                     /db_xref="GeneID:675"
     CDS             1..140
                     /gene="BRCA2"
                     /product="BRCA2 protein"
ORIGIN
        1 atgcctattg gatccaaaga gaggccaaca ttttttgaaa tttttaagac acgctgcaac
       61 aaagcagatt taggaccaat aagtcttaat tggtttgaag aactttcttc agaagctcca
      121 gggatcagaa tttggcacag
//
"""
    seqs = read_genbank(gb_content)
    assert len(seqs) == 1
    s = seqs[0]
    assert s.seq.upper().startswith("ATGCCTATTG")
    assert s.molecule_type == MoleculeType.DNA
    gene_features = [a for a in s.annotations if a.feature_type == "gene"]
    assert len(gene_features) >= 1
    assert gene_features[0].qualifiers.get("gene") == "BRCA2"


def test_write_fasta_line_wrap() -> None:
    from ..models.sequence import Sequence
    long_seq = "A" * 150
    seqs = [Sequence(id="test", seq=long_seq, molecule_type=MoleculeType.DNA)]
    out = write_fasta(seqs)
    lines = [ln for ln in out.splitlines() if not ln.startswith(">")]
    assert all(len(ln) <= 60 for ln in lines)
    assert "".join(lines) == long_seq


def test_parse_vcf_real(sample_vcf: str) -> None:
    variants = parse_vcf(sample_vcf)
    assert len(variants) == 4
    # BRCA2 pathogenic variant
    brca2 = next(v for v in variants if "BRCA2" in v.info.get("GENEINFO", ""))
    assert brca2.chrom == "13"
    assert brca2.pos == 32914437
    assert brca2.ref == "A"
    assert brca2.alt == ["C"]
    assert brca2.gene == "BRCA2"


def test_parse_vcf_types(sample_vcf: str) -> None:
    from ..models.variant import VariantType
    variants = parse_vcf(sample_vcf)
    # CFTR deletion (CTT→C)
    cftr = next(v for v in variants if "CFTR" in v.info.get("GENEINFO", ""))
    assert cftr.variant_type == VariantType.DEL


def test_write_vcf_roundtrip(sample_vcf: str) -> None:
    variants = parse_vcf(sample_vcf)
    vcf_out = write_vcf(variants)
    variants2 = parse_vcf(vcf_out)
    assert len(variants2) == len(variants)
    assert variants2[0].chrom == variants[0].chrom
    assert variants2[0].pos == variants[0].pos


def test_parse_gff3_real(sample_gff3: str) -> None:
    annotations = parse_gff3(sample_gff3)
    assert len(annotations) == 4
    gene = next(a for a in annotations if a.feature_type == "gene")
    assert gene.qualifiers.get("Name") == "BRCA1"
    assert gene.start == 43044294  # 0-based (44295-1)
    assert gene.end == 43125483
    from ..models.sequence import Strand
    assert gene.strand == Strand.MINUS


def test_gff3_roundtrip(sample_gff3: str) -> None:
    annotations = parse_gff3(sample_gff3)
    gff_out = write_gff3(annotations, seqid="chr17")
    annotations2 = parse_gff3(gff_out)
    assert len(annotations2) == len(annotations)
    types1 = {a.feature_type for a in annotations}
    types2 = {a.feature_type for a in annotations2}
    assert types1 == types2
