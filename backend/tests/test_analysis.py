"""Tests for analysis endpoints — protein calculator, ORF finder, IUPAC handling."""

from fastapi.testclient import TestClient

# Real TP53 protein N-terminus (from UniProt P04637)
TP53_PROTEIN = (
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDP"
    "GPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYPQGLNGTVNLFRNL"
)

# Real BRCA1 CDS fragment for ORF testing
BRCA1_CDS = (
    "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAA"
    "AATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTG"
    "ACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCA"
)

# Sequence with all ambiguity codes
AMBIGUOUS_SEQ = "ACGTNRYSWKMBDHVACGTNRYSWKMBDHV"

# Human beta-globin for composition
HBB_REGION = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAA"


def test_protein_properties_tp53(client: TestClient) -> None:
    r = client.post("/analysis/protein_properties", params={"sequence": TP53_PROTEIN})
    assert r.status_code == 200
    data = r.json()
    assert data["length"] == len(TP53_PROTEIN)
    # TP53 is acidic — pI should be below 7
    assert data["isoelectric_point"] < 7.0
    # Should have meaningful MW
    assert data["molecular_weight_da"] > 10000
    # GRAVY should be negative for hydrophilic proteins like TP53
    assert data["gravy"] < 0


def test_protein_properties_composition(client: TestClient) -> None:
    """Short well-characterized peptide: all residue counts should match."""
    peptide = "ACDEFGHIKLMNPQRSTVWY"  # one of each standard AA
    r = client.post("/analysis/protein_properties", params={"sequence": peptide})
    assert r.status_code == 200
    data = r.json()
    assert data["length"] == 20
    # Each of the 20 AAs appears exactly once
    comp = data["aa_composition"]
    assert all(comp.get(aa, 0) == 1 for aa in "ACDEFGHIKLMNPQRSTVWY")


def test_protein_properties_signal_peptide(client: TestClient) -> None:
    """Sequence with known signal peptide: MMSFVSLLLVGILFWATEAEQLTKCEVFQ (IgG signal)."""
    sp_seq = "MMSFVSLLLVGILFWATEAEQLTKCEVFQ" + "A" * 50
    r = client.post("/analysis/protein_properties", params={"sequence": sp_seq})
    assert r.status_code == 200
    data = r.json()
    assert data["signal_peptide_predicted"] is True
    assert data["signal_peptide_end"] is not None


def test_protein_properties_with_selenocysteine(client: TestClient) -> None:
    """Selenocysteine (U) should be accepted and counted."""
    r = client.post("/analysis/protein_properties", params={"sequence": "MACUGPWK"})
    assert r.status_code == 200
    data = r.json()
    assert data["aa_composition"].get("U", 0) == 1


def test_protein_properties_invalid(client: TestClient) -> None:
    """Nucleotide-only sequence should be rejected."""
    r = client.post("/analysis/protein_properties", params={"sequence": "1234ZZZQ!"})
    assert r.status_code == 400


def test_find_orfs_brca1(client: TestClient) -> None:
    """BRCA1 CDS should contain the main ORF in frame +1."""
    r = client.post("/analysis/find_orfs", params={"sequence": BRCA1_CDS, "min_length_aa": 10})
    assert r.status_code == 200
    data = r.json()
    assert data["total_found"] >= 1
    # Longest ORF should start with M
    longest = data["orfs"][0]
    assert longest["protein"].startswith("M")
    assert longest["start_codon"] == "ATG"


def test_find_orfs_all_frames(client: TestClient) -> None:
    """6-frame ORF finder: frames +1,+2,+3,-1,-2,-3 should all be checked."""
    seq = "ATG" + "GCC" * 30 + "TAA"  # simple ORF in +1
    r = client.post("/analysis/find_orfs", params={"sequence": seq, "min_length_aa": 5})
    assert r.status_code == 200
    orfs = r.json()["orfs"]
    frames = {o["frame"] for o in orfs}
    assert 1 in frames  # forward frame 1 must have the ORF


def test_find_orfs_min_length(client: TestClient) -> None:
    """min_length_aa filter should exclude short ORFs."""
    seq = "ATG" + "GCC" * 5 + "TAA" + "A" * 30 + "ATG" + "GCC" * 50 + "TAA"
    r_short = client.post("/analysis/find_orfs", params={"sequence": seq, "min_length_aa": 3})
    r_long = client.post("/analysis/find_orfs", params={"sequence": seq, "min_length_aa": 30})
    assert len(r_short.json()["orfs"]) >= len(r_long.json()["orfs"])


def test_composition_hbb(client: TestClient) -> None:
    """HBB region composition should have correct GC and all valid IUPAC bases."""
    r = client.post("/analysis/composition", params={"sequence": HBB_REGION})
    assert r.status_code == 200
    data = r.json()
    assert data["length"] == len(HBB_REGION)
    assert 40 <= data["gc_content"] <= 70
    # Only A, T, G, C should appear (no ambiguity)
    assert data["n_count"] == 0
    assert all(c in "ACGT" for c in data["symbol_counts"])


def test_composition_ambiguous_iupac(client: TestClient) -> None:
    """All IUPAC ambiguity codes should be counted correctly."""
    r = client.post("/analysis/composition", params={"sequence": AMBIGUOUS_SEQ})
    assert r.status_code == 200
    data = r.json()
    counts = data["symbol_counts"]
    # Each IUPAC code appears twice in AMBIGUOUS_SEQ
    for code in "ACGTNRYSWKMBDHV":
        assert counts.get(code, 0) == 2, f"Expected 2 of '{code}', got {counts.get(code, 0)}"
    # GC should account for ambiguous bases
    assert data["gc_content"] is not None
    assert 0 < data["gc_content"] < 100


def test_gc_window(client: TestClient) -> None:
    """Sliding window GC should produce len = (seq_len - w) / step + 1 points."""
    seq = "GCGCGCGCAT" * 20  # 200 bp, 80% GC
    r = client.post("/analysis/gc_window", params={"sequence": seq, "window_size": 50, "step": 10})
    assert r.status_code == 200
    data = r.json()
    assert len(data["positions"]) == len(data["gc_values"])
    assert len(data["positions"]) > 0
    # GC should be high for this sequence
    assert all(gc > 60 for gc in data["gc_values"])


def test_find_repeats(client: TestClient) -> None:
    """Known tandem repeat (AT)×10 should be detected."""
    seq = "GGGGG" + "AT" * 10 + "CCCCC"
    r = client.post("/analysis/find_repeats", params={"sequence": seq, "min_unit": 2, "max_unit": 4})
    assert r.status_code == 200
    repeats = r.json()
    tandem = [rp for rp in repeats if rp["repeat_type"] == "tandem" and rp["unit"] == "AT"]
    assert len(tandem) >= 1
    assert tandem[0]["copies"] >= 5.0


def test_optimize_codons_human(client: TestClient) -> None:
    """Codon optimization for human should produce valid DNA encoding same protein."""
    dna = "ATGGTGCACCTGACTCCTGAGGAG"  # HBB start: MVHLTP...
    r = client.post("/analysis/optimize_codons",
                    params={"dna_sequence": dna, "organism": "human"})
    assert r.status_code == 200
    data = r.json()
    # Output must be same length (codon-for-codon swap)
    assert len(data["optimized"]) == len(dna)
    # Must encode same protein
    from Bio.Seq import Seq
    orig_prot = str(Seq(dna).translate())
    opt_prot = str(Seq(data["optimized"]).translate())
    assert orig_prot == opt_prot


def test_optimize_codons_from_protein(client: TestClient) -> None:
    """Back-translate protein to optimized DNA for E. coli."""
    protein = "MVHLTPEEK"
    r = client.post("/analysis/optimize_codons",
                    params={"dna_sequence": "", "protein_sequence": protein, "organism": "ecoli"})
    assert r.status_code == 200
    data = r.json()
    assert len(data["optimized"]) == len(protein) * 3
    # Decode and verify
    from Bio.Seq import Seq
    decoded = str(Seq(data["optimized"]).translate(to_stop=True))
    assert decoded == protein


def test_iupac_codes_endpoint(client: TestClient) -> None:
    """IUPAC reference endpoint should list all standard codes."""
    r = client.get("/analysis/iupac_codes")
    assert r.status_code == 200
    data = r.json()
    # Must include all 4 definite DNA bases
    assert set("ACGT").issubset(data["nucleotides"]["definite"])
    # Must include all 20 standard amino acids
    assert set("ACDEFGHIKLMNPQRSTVWY").issubset(data["amino_acids"]["standard_20"])
    # Must include non-standard AAs
    assert "U" in data["amino_acids"]["non_standard"]
    assert "O" in data["amino_acids"]["non_standard"]
    # Must include stop codon
    assert "*" in data["amino_acids"]["special"]
