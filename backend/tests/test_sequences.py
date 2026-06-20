"""Tests for sequence editing endpoints using real gene sequences."""

from fastapi.testclient import TestClient


def test_add_and_get_sequence(client: TestClient) -> None:
    seq_data = {
        "id": "HBB_test",
        "name": "HBB",
        "description": "Human beta-globin mRNA (test)",
        "seq": "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAA",
        "molecule_type": "DNA",
    }
    r = client.post("/sequences/", json=seq_data)
    assert r.status_code == 201
    data = r.json()
    assert data["id"] == "HBB_test"
    assert data["length"] == 69

    r2 = client.get("/sequences/HBB_test")
    assert r2.status_code == 200
    assert r2.json()["seq"] == seq_data["seq"]


def test_list_sequences(client: TestClient) -> None:
    client.post("/sequences/", json={"id": "seq_list_test", "seq": "ATCG", "molecule_type": "DNA"})
    r = client.get("/sequences/")
    assert r.status_code == 200
    assert any(s["id"] == "seq_list_test" for s in r.json())


def test_delete_sequence(client: TestClient) -> None:
    client.post("/sequences/", json={"id": "to_delete", "seq": "ATCG", "molecule_type": "DNA"})
    r = client.delete("/sequences/to_delete")
    assert r.status_code == 204
    r2 = client.get("/sequences/to_delete")
    assert r2.status_code == 404


def test_reverse_complement(client: TestClient) -> None:
    # BRCA2 codon: ATG → CAT on minus strand
    client.post("/sequences/", json={"id": "rc_test", "seq": "ATGCCTATTGGATCC", "molecule_type": "DNA"})
    r = client.post("/sequences/rc_test/edit", json={"operation": "reverse_complement"})
    assert r.status_code == 200
    data = r.json()
    assert data["result_seq"] == "GGATCCAATAGGCAT"


def test_translate(client: TestClient) -> None:
    # ATG GTG CAC = Met Val His
    client.post("/sequences/", json={"id": "translate_test", "seq": "ATGGTGCAC", "molecule_type": "DNA"})
    r = client.post("/sequences/translate_test/edit", json={"operation": "translate"})
    assert r.status_code == 200
    assert r.json()["result_seq"] == "MVH"


def test_translate_hbb(client: TestClient, hbb_fasta: str) -> None:
    """Translate real HBB CDS — must match actual beta-globin N-terminus."""
    hbb_cds = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAA"
    client.post("/sequences/", json={"id": "hbb_cds2", "seq": hbb_cds, "molecule_type": "DNA"})
    r = client.post("/sequences/hbb_cds2/edit", json={"operation": "translate"})
    assert r.status_code == 200
    protein = r.json()["result_seq"]
    # Real HBB N-terminus: Met-Val-His-Leu-Thr-Pro-Glu-Glu-Lys-Ser-Ala-Val-Thr-Ala-Leu-Trp-Gly-Lys-Val-Asn-Val-Asp-Glu
    assert protein == "MVHLTPEEKSAVTALWGKVNVDE"


def test_translate_brca2(client: TestClient, brca2_fasta: str) -> None:
    """Translate real BRCA2 partial CDS — must match known N-terminus."""
    from backend.formats import read_fasta
    seqs = read_fasta(brca2_fasta)
    brca2_seq = seqs[0].seq
    client.post("/sequences/", json={"id": "brca2_t", "seq": brca2_seq, "molecule_type": "DNA"})
    r = client.post("/sequences/brca2_t/edit", json={"operation": "translate"})
    assert r.status_code == 200
    protein = r.json()["result_seq"]
    # BRCA2 N-terminus (from NM_000059.4 partial): MPIGSKERPTFFEIFKTRCNKADLGPISLNWFEELSSEAPGIRIWH
    assert protein.startswith("MPIGS"), f"Expected MPIGS…, got {protein[:10]}…"
    assert "KERCNKADLG" not in protein  # sanity: not a garbled sequence


def test_insert(client: TestClient) -> None:
    client.post("/sequences/", json={"id": "ins_test", "seq": "AAATTTGGG", "molecule_type": "DNA"})
    r = client.post("/sequences/ins_test/edit", json={
        "operation": "insert", "position": 3, "insert_seq": "CCC"
    })
    assert r.status_code == 200
    assert r.json()["result_seq"] == "AAACCCTTTGGG"


def test_delete(client: TestClient) -> None:
    client.post("/sequences/", json={"id": "del_test", "seq": "AAATTTGGG", "molecule_type": "DNA"})
    r = client.post("/sequences/del_test/edit", json={
        "operation": "delete", "position": 3, "end_position": 6
    })
    assert r.status_code == 200
    assert r.json()["result_seq"] == "AAAGGG"


def test_replace(client: TestClient) -> None:
    client.post("/sequences/", json={"id": "rep_test", "seq": "AAATTTGGG", "molecule_type": "DNA"})
    r = client.post("/sequences/rep_test/edit", json={
        "operation": "replace", "position": 3, "end_position": 6, "replacement": "CCC"
    })
    assert r.status_code == 200
    assert r.json()["result_seq"] == "AAACCCGGG"


def test_gc_content(client: TestClient) -> None:
    client.post("/sequences/", json={"id": "gc_test", "seq": "ATGCATGC", "molecule_type": "DNA"})
    r = client.get("/sequences/gc_test/gc_content")
    assert r.status_code == 200
    data = r.json()
    assert data["gc_content"] == 50.0
    assert data["length"] == 8


def test_gc_content_hbb(client: TestClient) -> None:
    """Real HBB sequence GC content should be within biologically expected range."""
    hbb = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAA"
    client.post("/sequences/", json={"id": "hbb_gc", "seq": hbb, "molecule_type": "DNA"})
    r = client.get("/sequences/hbb_gc/gc_content")
    assert r.status_code == 200
    gc = r.json()["gc_content"]
    assert 40.0 <= gc <= 70.0  # typical coding sequence


def test_codon_usage(client: TestClient) -> None:
    client.post("/sequences/", json={"id": "codon_test", "seq": "ATGATGATG", "molecule_type": "DNA"})
    r = client.get("/sequences/codon_test/codon_usage")
    assert r.status_code == 200
    data = r.json()
    assert data["total_codons"] == 3
    assert data["codon_usage"]["ATG"]["count"] == 3


def test_find_motif_kozak(client: TestClient) -> None:
    """Find Kozak consensus (GCCACCATG) in HBB — real start codon context."""
    seq = "GCCACCATGGTGCACCTGACTCCT"
    client.post("/sequences/", json={"id": "kozak_test", "seq": seq, "molecule_type": "DNA"})
    r = client.get("/sequences/kozak_test/find_motif?motif=ATG")
    assert r.status_code == 200
    data = r.json()
    assert data["count"] == 1
    assert 6 in data["positions"]


def test_find_motif_iupac(client: TestClient) -> None:
    """IUPAC N matches any base."""
    seq = "AACGTATCGATCG"
    client.post("/sequences/", json={"id": "iupac_test", "seq": seq, "molecule_type": "DNA"})
    r = client.get("/sequences/iupac_test/find_motif?motif=NNN")
    assert r.status_code == 200
    assert r.json()["count"] > 0


def test_invalid_operation(client: TestClient) -> None:
    client.post("/sequences/", json={"id": "op_test", "seq": "ATCG", "molecule_type": "DNA"})
    r = client.post("/sequences/op_test/edit", json={"operation": "nonsense"})
    assert r.status_code == 400


def test_transcribe(client: TestClient) -> None:
    client.post("/sequences/", json={"id": "tx_test", "seq": "ATCGATCG", "molecule_type": "DNA"})
    r = client.post("/sequences/tx_test/edit", json={"operation": "transcribe"})
    assert r.status_code == 200
    assert r.json()["result_seq"] == "AUCGAUCG"


def test_complement(client: TestClient) -> None:
    client.post("/sequences/", json={"id": "comp_test", "seq": "ATCG", "molecule_type": "DNA"})
    r = client.post("/sequences/comp_test/edit", json={"operation": "complement"})
    assert r.status_code == 200
    assert r.json()["result_seq"] == "TAGC"
