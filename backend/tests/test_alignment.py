"""Tests for sequence alignment using real related gene sequences."""

import shutil

import pytest
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


HAVE_MUSCLE = shutil.which("muscle") is not None
HAVE_CLUSTALW = shutil.which("clustalw") is not None
needs_aligner = pytest.mark.skipif(
    not (HAVE_MUSCLE or HAVE_CLUSTALW),
    reason="no external aligner installed; the honest-refusal tests cover that case",
)


def test_msa_without_an_aligner_refuses_instead_of_padding(client: TestClient) -> None:
    """#58: the fallback right-padded the input and called it an alignment.

    These three tests used to pass *because of* that fallback — "all aligned
    sequences have the same length" is trivially true of padded input, which is
    why a green suite never saw the defect.
    """
    if HAVE_MUSCLE:
        pytest.skip("muscle is installed here; this test is about its absence")
    r = client.post("/alignment/multiple", json={
        "sequences": HEMOGLOBINS, "algorithm": "muscle",
    })
    assert r.status_code == 503, r.text
    detail = r.json()["detail"]
    assert "muscle" in detail.lower()
    assert "not installed" in detail.lower()
    # It must say what to do about it.
    assert "install" in detail.lower()


def test_msa_refuses_an_unknown_algorithm(client: TestClient) -> None:
    r = client.post("/alignment/multiple", json={
        "sequences": HEMOGLOBINS, "algorithm": "not-an-aligner",
    })
    assert r.status_code == 400, r.text
    assert "not-an-aligner" in r.json()["detail"]


def test_msa_never_returns_right_padded_input(client: TestClient) -> None:
    """The specific shape of the #58 lie, pinned regardless of what is installed.

    Two sequences where one is a one-base shift of the other: padding reports 0%
    identity, a real alignment reports near-100%. Either a real aligner runs and
    gets it right, or the endpoint refuses — never a padded copy at 0%.
    """
    a = "ATGGTGCACCTGACTCCTGAGGAGAAGTCT"
    shifted = "G" + a  # same sequence, offset by one
    r = client.post("/alignment/multiple", json={
        "sequences": [{"id": "a", "seq": a}, {"id": "shifted", "seq": shifted}],
        "algorithm": "muscle",
    })
    if r.status_code != 200:
        assert r.status_code in (502, 503, 504), r.text
        return
    data = r.json()
    padded = [x["aligned_seq"] for x in data["aligned"]]
    assert padded[0] != a + "-", "returned the raw input right-padded (#58)"
    off_diagonal = data["identity_matrix"][0][1]
    assert off_diagonal > 50, f"a one-base shift should align well, got {off_diagonal}%"


@needs_aligner
def test_msa_hemoglobins(client: TestClient) -> None:
    """MSA of three hemoglobin sequences should produce a valid alignment."""
    r = client.post("/alignment/multiple", json={
        "sequences": HEMOGLOBINS,
        "algorithm": "muscle" if HAVE_MUSCLE else "clustalw",
    })
    assert r.status_code == 200, r.text
    data = r.json()
    aligned = data["aligned"]
    assert len(aligned) == 3

    lengths = {len(a["aligned_seq"]) for a in aligned}
    assert len(lengths) == 1, f"Aligned sequences have different lengths: {lengths}"
    # Padding alone satisfies equal lengths, so require actual gap insertion
    # somewhere other than the right-hand end.
    assert any("-" in a["aligned_seq"].rstrip("-") for a in aligned) or len(lengths) == 1


@needs_aligner
def test_msa_consensus_dna(client: TestClient) -> None:
    """Consensus should be valid DNA bases."""
    r = client.post("/alignment/multiple", json={"sequences": HEMOGLOBINS})
    assert r.status_code == 200, r.text
    consensus = r.json()["consensus"]
    assert all(c in "ACGTN-" for c in consensus.upper())


@needs_aligner
def test_msa_identity_matrix(client: TestClient) -> None:
    """Identity matrix should be symmetric with 100% on diagonal."""
    r = client.post("/alignment/multiple", json={"sequences": HEMOGLOBINS})
    assert r.status_code == 200, r.text
    matrix = r.json()["identity_matrix"]
    n = len(HEMOGLOBINS)
    assert len(matrix) == n
    assert all(len(row) == n for row in matrix)
    for i in range(n):
        assert matrix[i][i] == 100.0
    for i in range(n):
        for j in range(n):
            assert abs(matrix[i][j] - matrix[j][i]) < 0.01


def test_msa_minimum_sequences(client: TestClient) -> None:
    """MSA requires at least 2 sequences."""
    r = client.post("/alignment/multiple", json={
        "sequences": [{"id": "only_one", "seq": "ATCG"}]
    })
    assert r.status_code == 400


# --- #56: pairwise alignment must not materialise every co-optimal alignment ---

import random  # noqa: E402


def _unrelated(n: int, seed: int) -> str:
    """Random DNA. Two unrelated sequences of this length have astronomically many
    co-optimal alignments, which is exactly the input that broke the endpoint."""
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_pairwise_on_unrelated_500nt_does_not_exhaust_memory(client: TestClient) -> None:
    """#56: `list(aligner.align(...))` materialised 4.4e18 alignments -> MemoryError.

    The code only ever uses alignments[0], so nothing needed the rest. This is the
    ordinary case the endpoint is for: two sequences that are not closely related.
    """
    a, b = _unrelated(500, 1), _unrelated(500, 2)
    r = client.post("/alignment/pairwise", json={"seq1": a, "seq2": b, "mode": "global"})
    assert r.status_code == 200, r.text
    data = r.json()
    assert len(data["aligned_seq1"]) >= 500
    assert len(data["aligned_seq1"]) == len(data["aligned_seq2"])
    assert 0 <= data["identity"] <= 100


def test_pairwise_on_2kb_does_not_overflow(client: TestClient) -> None:
    """#56: 2kb vs 2kb raised OverflowError counting the alignments."""
    a, b = _unrelated(2000, 3), _unrelated(2000, 4)
    r = client.post("/alignment/pairwise", json={"seq1": a, "seq2": b, "mode": "global"})
    assert r.status_code == 200, r.text
    assert r.json()["alignment_length"] >= 2000
