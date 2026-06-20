"""Tests for sequence alignment using real related gene sequences."""

from fastapi.testclient import TestClient

# Real hemoglobin alpha vs beta (related but distinct)
HBA1 = "ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGG"
HBB = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGG"

# Human vs chimp TP53 exon (highly similar)
TP53_HUMAN = "ATGGAGGAGCCGCAGTCAGATCCTAGCGTTCGAGTCCTGCCTACTGGCCTGCACC"
TP53_CHIMP = "ATGGAGGAGCCGCAGTCAGATCCTAGCGTTCAAGTCCTGCCTACTGGCCTGCACC"  # 1 SNP

# Multiple hemoglobin sequences for MSA
HEMOGLOBINS = [
    {"id": "HBB_human", "seq": "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCC"},
    {"id": "HBA1_human", "seq": "ATGGTGCTGTCTCCTGCCGACAAGACCAACGTC"},
    {"id": "HBB_mouse", "seq": "ATGGTGCACCTGACTCCTGAAGAGAAGGCTGCC"},
]


def test_pairwise_global_alignment(client: TestClient) -> None:
    r = client.post("/alignment/pairwise", json={
        "seq1": HBA1,
        "seq2": HBB,
        "mode": "global",
    })
    assert r.status_code == 200
    data = r.json()
    assert data["score"] > 0
    assert len(data["aligned_seq1"]) > 0
    assert len(data["aligned_seq1"]) == len(data["aligned_seq2"])
    assert 0 < data["identity"] <= 100
    assert 0 < data["similarity"] <= 100
    assert data["alignment_length"] > 0


def test_pairwise_near_identical(client: TestClient) -> None:
    """Human vs chimp TP53 — should have very high identity (differ by 1 nt)."""
    r = client.post("/alignment/pairwise", json={
        "seq1": TP53_HUMAN,
        "seq2": TP53_CHIMP,
        "mode": "global",
    })
    assert r.status_code == 200
    data = r.json()
    # 1 SNP in ~56 nt → identity >98%
    assert data["identity"] > 95.0


def test_pairwise_local_alignment(client: TestClient) -> None:
    """Local alignment should find a high-scoring sub-region."""
    r = client.post("/alignment/pairwise", json={
        "seq1": "AAAAATGGTGCACATTT",
        "seq2": "GGGGATGGTGCACGGG",
        "mode": "local",
    })
    assert r.status_code == 200
    data = r.json()
    assert data["score"] > 0


def test_pairwise_self_alignment(client: TestClient) -> None:
    """Self-alignment should yield 100% identity."""
    seq = "ATGCGATCGATCG"
    r = client.post("/alignment/pairwise", json={"seq1": seq, "seq2": seq, "mode": "global"})
    assert r.status_code == 200
    assert r.json()["identity"] == 100.0


def test_msa_hemoglobins(client: TestClient) -> None:
    """MSA of three hemoglobin sequences should produce valid alignment."""
    r = client.post("/alignment/multiple", json={
        "sequences": HEMOGLOBINS,
        "algorithm": "muscle",
    })
    assert r.status_code == 200
    data = r.json()
    aligned = data["aligned"]
    assert len(aligned) == 3

    # All aligned sequences should have the same length
    lengths = {len(a["aligned_seq"]) for a in aligned}
    assert len(lengths) == 1, f"Aligned sequences have different lengths: {lengths}"


def test_msa_consensus_dna(client: TestClient) -> None:
    """Consensus should be valid DNA bases."""
    r = client.post("/alignment/multiple", json={"sequences": HEMOGLOBINS})
    assert r.status_code == 200
    consensus = r.json()["consensus"]
    assert all(c in "ACGTN-" for c in consensus.upper())


def test_msa_identity_matrix(client: TestClient) -> None:
    """Identity matrix should be symmetric with 100% on diagonal."""
    r = client.post("/alignment/multiple", json={"sequences": HEMOGLOBINS})
    assert r.status_code == 200
    matrix = r.json()["identity_matrix"]
    n = len(HEMOGLOBINS)
    assert len(matrix) == n
    assert all(len(row) == n for row in matrix)
    # Diagonal = 100
    for i in range(n):
        assert matrix[i][i] == 100.0
    # Symmetry
    for i in range(n):
        for j in range(n):
            assert abs(matrix[i][j] - matrix[j][i]) < 0.01


def test_msa_minimum_sequences(client: TestClient) -> None:
    """MSA requires at least 2 sequences."""
    r = client.post("/alignment/multiple", json={
        "sequences": [{"id": "only_one", "seq": "ATCG"}]
    })
    assert r.status_code == 400
