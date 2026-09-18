#!/usr/bin/env python3
"""
Phase 1+2: Discovery & Analysis
Fetches real biological data from NCBI/Ensembl and builds a test corpus
with expected output shapes for every Oligolia endpoint.

Output: .claude/qa/corpus/corpus.json
"""

import json
import time
import urllib.request
import urllib.parse
import urllib.error
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import _state  # noqa: E402

BASE = Path(__file__).parent
CORPUS_DIR = BASE / "corpus"
CORPUS_DIR.mkdir(exist_ok=True)

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL_BASE = "https://rest.ensembl.org"


def get(url, headers=None, retries=3, delay=0.4):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers or {"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=15) as r:
                return json.loads(r.read())
        except Exception as e:
            if i == retries - 1:
                print(f"  WARN: {url} failed: {e}")
                return None
            time.sleep(delay * (i + 1))


def post(url, payload, headers=None):
    data = json.dumps(payload).encode()
    hdrs = {"Content-Type": "application/json", **(headers or {})}
    try:
        req = urllib.request.Request(url, data=data, headers=hdrs, method="POST")
        with urllib.request.urlopen(req, timeout=15) as r:
            return json.loads(r.read())
    except Exception as e:
        print(f"  WARN: POST {url} failed: {e}")
        return None


# ── Real gene sequences from NCBI ────────────────────────────────────────────

GENE_QUERIES = [
    # (gene_symbol, ncbi_gene_id, description)
    ("BRCA1", "672",  "Breast cancer susceptibility — well-known CRISPR target"),
    ("TP53",  "7157", "Tumour suppressor — high variant density"),
    ("HBB",   "3043", "Haemoglobin beta — classic sickle-cell gene"),
    ("CFTR",  "1080", "Cystic fibrosis — large gene, many exons"),
    ("EGFR",  "1956", "Oncogene — common therapy target"),
]

def fetch_ncbi_sequence(gene_id: str, gene_symbol: str) -> dict | None:
    """Fetch NM_ RefSeq mRNA sequence for a gene via NCBI eutils."""
    print(f"  Fetching {gene_symbol} (gene_id={gene_id}) from NCBI...")
    # Search for RefSeq mRNA accession
    search_url = (
        f"{NCBI_BASE}/esearch.fcgi?db=nucleotide&term={gene_symbol}[Gene+Name]+AND+"
        f"Homo+sapiens[Organism]+AND+mRNA[Filter]+AND+RefSeq[Filter]"
        f"&retmax=1&retmode=json"
    )
    result = get(search_url)
    if not result:
        return None
    ids = result.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None

    # Fetch sequence (first 2000 bp for speed)
    fetch_url = (
        f"{NCBI_BASE}/efetch.fcgi?db=nucleotide&id={ids[0]}"
        f"&rettype=fasta&retmode=text&seq_start=1&seq_stop=2000"
    )
    try:
        req = urllib.request.Request(fetch_url)
        with urllib.request.urlopen(req, timeout=15) as r:
            fasta = r.read().decode("utf-8")
        lines = fasta.strip().split("\n")
        header = lines[0]
        seq = "".join(lines[1:]).upper()
        return {"header": header, "seq": seq, "gene_symbol": gene_symbol, "ncbi_id": ids[0]}
    except Exception as e:
        print(f"  WARN: efetch for {gene_symbol} failed: {e}")
        return None
    finally:
        time.sleep(0.35)  # NCBI rate limit: 3 req/s unauthenticated


def fetch_ensembl_sequence(gene_symbol: str) -> dict | None:
    """Fetch genomic sequence from Ensembl REST."""
    print(f"  Fetching {gene_symbol} from Ensembl...")
    # Lookup gene ID
    lookup = get(
        f"{ENSEMBL_BASE}/xrefs/symbol/homo_sapiens/{gene_symbol}?content-type=application/json"
    )
    if not lookup:
        return None
    gene_ids = [x["id"] for x in lookup if x.get("id", "").startswith("ENSG")]
    if not gene_ids:
        return None

    ensg = gene_ids[0]
    # Fetch sequence (first 1500 bp)
    seq_data = get(
        f"{ENSEMBL_BASE}/sequence/id/{ensg}?content-type=application/json&type=genomic&expand_5prime=0&expand_3prime=0"
    )
    if not seq_data or "seq" not in seq_data:
        return None
    seq = seq_data["seq"][:1500].upper()
    return {"ensembl_id": ensg, "gene_symbol": gene_symbol, "seq": seq}


# ── Edge case sequences (no network needed) ───────────────────────────────────

EDGE_CASES = [
    {
        "id": "edge_short",
        "desc": "Very short sequence — 8 bp, below primer design minimum",
        "seq": "ATGCATGC",
        "expect_primer_fail": True,
        "expect_crispr_fail": True,
    },
    {
        "id": "edge_ambiguous",
        "desc": "Sequence with IUPAC ambiguity codes (R, Y, S, W, K, M)",
        "seq": "ATGRYGCATGCSWKMATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC",
        "expect_crispr_partial": True,
    },
    {
        "id": "edge_lowercase",
        "desc": "Lowercase sequence — should be normalised",
        "seq": "atggtgcacctgactcctgaggagaagtctgccgttactgccctgtggggcaaggtgaacgtg",
        "expect_normalised": True,
    },
    {
        "id": "edge_rna",
        "desc": "RNA sequence (U instead of T)",
        "seq": "AUGGUGCACCUGACUCCUGAGGAGAAGUCUGCCGUUACUGCCCUGUGGGGCAAGGUGAACGUG",
        "seq_type": "RNA",
    },
    {
        "id": "edge_protein",
        "desc": "Protein sequence — single-letter AA codes",
        "seq": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
        "seq_type": "protein",
    },
    {
        "id": "edge_multiline_fasta",
        "desc": "Multi-sequence FASTA — for MSA endpoint",
        "fasta": ">HBB_human\nATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTG\n>HBB_chimp\nATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTG\n>HBB_mouse\nATGGTGCACCTGACCCCTGAGGAGAAGTCCGCCGTTACTGCCCTGTGGGGCAAGGTG",
    },
    {
        "id": "edge_vcf_multialt",
        "desc": "VCF with multi-ALT allele — known annotation bug",
        "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t925952\t.\tG\tA,T\t.\tPASS\t.\nchr17\t43092919\t.\tA\tG\t50\tPASS\t.\n",
    },
    {
        "id": "edge_vcf_no_chr_prefix",
        "desc": "VCF using bare chromosome numbers (1, not chr1) — common from GATK",
        "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t925952\t.\tG\tA\t.\tPASS\t.\n17\t43092919\t.\tA\tG\t50\tPASS\t.\n",
    },
]


# ── Real variant data ─────────────────────────────────────────────────────────

REAL_VARIANTS = [
    # (chrom, pos, ref, alt, description)
    ("chr17", "43092919", "A",  "G",   "BRCA1 known pathogenic — rs28897672"),
    ("chr11", "5246945",  "T",  "A",   "HBB sickle-cell — rs334"),
    ("chr7",  "117559590","CTT","C",   "CFTR ΔF508 — most common CF variant"),
    ("chr17", "7674220",  "C",  "T",   "TP53 hotspot — rs28934578"),
    ("chr7",  "55174014", "C",  "T",   "EGFR L858R — common NSCLC driver"),
    ("chr1",  "925952",   "G",  "A,T", "Multi-ALT — expected to fail annotation (bug #1)"),
]


# ── MSA sequences (hemoglobin family for meaningful alignment) ────────────────

MSA_SEQS = [
    {"id": "HBA1_human",  "seq": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"},
    {"id": "HBB_human",   "seq": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"},
    {"id": "HBG1_human",  "seq": "MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGAFAQLSELHCDKLHVDPENFRLLGNVLVTVLAHHFGKEFTPELQASYQKVVAGVANALAHGYH"},
    {"id": "HBE1_human",  "seq": "MVHFTAEEKAAVTSLWSKMNVEEAGGEALGRLLVVYPWTQRFFDSFGNLSSPSAILGNPKVKAHGKKVLTSFGDAIKNMDNLKPAFAKLSELHCDKLHVDPENFRLLGNVLVCVLARNFGKEFTPQMQAAYQKVVAGVANALAHRYH"},
]


# ── Build corpus ──────────────────────────────────────────────────────────────

def build_corpus() -> dict:
    print("\n=== Phase 1: Discovery — fetching real biological data ===\n")

    gene_sequences = []
    skipped: list[dict] = []
    for symbol, ncbi_id, desc in GENE_QUERIES:
        seq_data = fetch_ncbi_sequence(ncbi_id, symbol)
        if seq_data:
            seq_data["description"] = desc
            gene_sequences.append(seq_data)
            print(f"  ✓ {symbol}: {len(seq_data['seq'])} bp")
        else:
            # NO SUBSTITUTE SEQUENCE. This used to append 92 nt of HBB exon 1
            # under the requested gene's symbol, so a failed BRCA1 fetch produced
            # a corpus entry labelled BRCA1 containing HBB -- and since `expected`
            # is derived from whatever sequence is present, every assertion passed
            # and the run reported green (#73).
            #
            # Scale mattered too: 92 nt instead of ~2000 bp, which is the whole
            # reason for pulling real genes.
            print(f"  ✗ {symbol}: fetch failed — SKIPPING (no substitute sequence)")
            skipped.append({"gene_symbol": symbol, "reason": "NCBI fetch failed"})

    if skipped:
        print(f"\n  {len(skipped)} gene(s) skipped: "
              + ", ".join(s["gene_symbol"] for s in skipped))
    if not gene_sequences:
        raise SystemExit(
            "scout: every NCBI fetch failed, so there is no corpus to build. "
            "Refusing to emit one rather than substituting sequences (#73)."
        )

    print("\n=== Phase 2: Analysis — defining expected output shapes ===\n")

    corpus = {
        "version": "1.0",
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "test_cases": [],
    }

    # ── Sequences CRUD ────────────────────────────────────────────────────────
    for gs in gene_sequences:
        corpus["test_cases"].append({
            "id": f"seq_create_{gs['gene_symbol']}",
            "phase": "execution",
            # Trailing slash: FastAPI 307-redirects POST /sequences, which the runner
            # recorded as a failure (5 of 22).
            "endpoint": "POST /sequences/",
            "payload": {
                "id": gs["gene_symbol"].lower(),
                "name": gs["gene_symbol"],
                "seq": gs["seq"],
                "seq_type": "DNA",
                "organism": "Homo sapiens",
                "description": gs.get("description", ""),
            },
            "expected": {
                "status": 201,
                "body_contains": ["id", "seq", "length"],
                "length": len(gs["seq"]),
            },
            "tags": ["sequences", "crud"],
        })

    # ── CRISPR design ─────────────────────────────────────────────────────────
    for gs in gene_sequences[:3]:  # BRCA1, TP53, HBB — well-validated targets
        # The CasType enum values, not lowercase nicknames. "cas9" produced a
        # 422 on every single CRISPR case -- 10 of the 22 failures in the first
        # real run, all of them corpus defects rather than app defects (#69).
        for cas in ["SpCas9", "SpCas9-HF1", "AsCas12a"]:
            corpus["test_cases"].append({
                "id": f"crispr_{gs['gene_symbol']}_{cas}",
                "phase": "execution",
                "endpoint": "POST /crispr/design",
                "payload": {
                    "target_sequence": gs["seq"][:500],  # first 500 bp
                    "cas_type": cas,
                    "guide_length": 20,
                    "max_guides": 5,
                    "check_off_targets": False,
                },
                "expected": {
                    "status": 200,
                    "body_contains": ["guides", "total_candidates"],
                    "guides_min": 0,  # some sequences may have no PAM
                    "each_guide_has": ["sequence", "position", "strand", "score"],
                    "score_range": [0.0, 1.0],
                },
                "tags": ["crispr", cas],
            })

    # ── CRISPR edge cases ─────────────────────────────────────────────────────
    corpus["test_cases"].append({
        "id": "crispr_edge_short",
        "endpoint": "POST /crispr/design",
        "payload": {"target_sequence": "ATGCATGC", "cas_type": "SpCas9", "guide_length": 20, "max_guides": 5},
        "expected": {"status": [400, 422], "body_contains": ["detail"]},
        "tags": ["crispr", "edge_case", "negative"],
    })
    corpus["test_cases"].append({
        # Was `crispr_edge_inverted_tm` with a 15 nt target and an expectation of
        # 200. It is about ambiguous bases, not Tm, and 15 nt is below the 23 nt
        # SpCas9 minimum -- so it asserted 200 against a correct 400 and reported
        # the guard as a defect. Long enough to actually exercise guide-finding
        # now, with an N-run in the middle where a guide would otherwise sit.
        "id": "crispr_edge_ambiguous_bases",
        "endpoint": "POST /crispr/design",
        "payload": {
            "target_sequence": "ATGGCCTGTGGGCATTTGGCCAANNNNNNNTTAGGCCATGGACGTGGCATCACGTGG",
            "cas_type": "SpCas9", "max_guides": 5,
        },
        "expected": {
            "status": 200,
            "body_contains": ["guides", "total_candidates"],
        },
        "tags": ["crispr", "edge_case", "ambiguous_bases"],
    })

    # ── Primer design ─────────────────────────────────────────────────────────
    for gs in gene_sequences:
        corpus["test_cases"].append({
            "id": f"primers_{gs['gene_symbol']}_standard",
            "endpoint": "POST /primers/design",
            "payload": {
                "template": gs["seq"],
                "product_min": 100,
                "product_max": 500,
                "primer_len_min": 18,
                "primer_len_max": 24,
                "tm_min": 55.0,
                "tm_max": 65.0,
                "gc_min": 40.0,
                "gc_max": 70.0,
                "max_pairs": 5,
            },
            "expected": {
                "status": 200,
                "is_list": True,
                # PrimerPair nests the primers: forward.tm / reverse.tm. There is no
                # tm_fwd or tm_rev and never was, so this asked for fields the API
                # does not have and reported their absence as a defect (#69).
                "each_item_has": ["forward", "reverse", "product_size", "penalty"],
                "product_size_range": [100, 500],
                "tm_range": [55.0, 65.0],
            },
            "tags": ["primers", "standard"],
        })

    # Inverted Tm range — should error or return empty with warning (currently silent)
    corpus["test_cases"].append({
        "id": "primers_edge_inverted_tm",
        "endpoint": "POST /primers/design",
        "payload": {
            "template": gene_sequences[0]["seq"],
            "tm_min": 65.0,
            "tm_max": 55.0,  # inverted
            "product_min": 100,
            "product_max": 500,
        },
        "expected": {
            "status": [200, 400, 422],
            "note": "Bug #2: Currently returns 200 with empty list, no error. Should be 400.",
            "known_bug": "primers_inverted_tm_silent_failure",
        },
        "tags": ["primers", "edge_case", "negative", "known_bug"],
    })

    # ── Pairwise alignment ────────────────────────────────────────────────────
    corpus["test_cases"].append({
        "id": "align_pairwise_nw_similar",
        "endpoint": "POST /alignment/pairwise",
        "payload": {
            "seq1": "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCC",
            "seq2": "ATGGTGCACCTGACCCCTGAGGAGAAGTCCGCC",
            "mode": "global",
            "match_score": 1,
            "mismatch_score": -1,
            "gap_open": -2,
            "gap_extend": -0.5,
        },
        "expected": {
            "status": 200,
            "body_contains": ["score", "identity", "aligned_seq1", "aligned_seq2"],
            "identity_range": [80.0, 100.0],
        },
        "tags": ["alignment", "pairwise"],
    })
    corpus["test_cases"].append({
        "id": "align_pairwise_sw_local",
        "endpoint": "POST /alignment/pairwise",
        "payload": {
            "seq1": gene_sequences[0]["seq"][:200],
            "seq2": gene_sequences[0]["seq"][100:300],
            "mode": "local",
        },
        "expected": {
            "status": 200,
            "body_contains": ["score", "identity"],
            "identity_range": [0.0, 100.0],
        },
        "tags": ["alignment", "pairwise", "local"],
    })

    # ── Multiple sequence alignment ───────────────────────────────────────────
    corpus["test_cases"].append({
        "id": "align_msa_hemoglobin",
        "endpoint": "POST /alignment/multiple",
        "payload": {"sequences": MSA_SEQS, "algorithm": "muscle"},
        "expected": {
            # 200 with an aligner installed; 503 without one, which is the honest
            # refusal from #58 rather than a padded fake alignment.
            "status": 200,
            "accept_status": [200, 503],
            "body_contains": ["aligned", "consensus", "identity_matrix"],
            "aligned_count": len(MSA_SEQS),
            "matrix_shape": [len(MSA_SEQS), len(MSA_SEQS)],
        },
        "tags": ["alignment", "msa", "protein"],
    })

    # ── Variant annotation ────────────────────────────────────────────────────
    variants_payload = [
        # `alt` is a list of alleles, not a string.
        {"chrom": c, "pos": int(p), "ref": r, "alt": [a], "id": f"{c}-{p}-{r}-{a}"}
        for c, p, r, a, _ in REAL_VARIANTS
    ]
    corpus["test_cases"].append({
        "id": "variants_annotate_real",
        "endpoint": "POST /variants/annotate",
        "payload": {
            "variants": variants_payload,
            "annotate_clinvar": True,
            "annotate_gnomad": True,
        },
        "expected": {
            "status": 200,
            "body_contains": ["variants", "total", "annotated"],
            "total": len(variants_payload),
            "note": "Multi-ALT variant (chr1-925952-G-A,T) expected to fail gnomAD lookup — bug #1",
        },
        "tags": ["variants", "annotation", "known_bug"],
    })

    # ── VCF parse ─────────────────────────────────────────────────────────────
    for ec in EDGE_CASES:
        if "vcf" in ec:
            corpus["test_cases"].append({
                "id": f"vcf_parse_{ec['id']}",
                "endpoint": "POST /files/parse/vcf",
                "payload_raw": ec["vcf"],
                "payload_type": "vcf_upload",
                "expected": {
                    "status": 200,
                    "is_list": True,
                    "note": ec["desc"],
                },
                "tags": ["vcf", "files", "edge_case"],
            })

    # ── Sequence analysis ─────────────────────────────────────────────────────
    for gs in gene_sequences[:2]:
        corpus["test_cases"].append({
            "id": f"analysis_orfs_{gs['gene_symbol']}",
            "endpoint": "POST /analysis/find_orfs",
            # `sequence` is a query parameter on this route, not a body field.
            "query": {"sequence": gs["seq"], "min_length_aa": 30},
            "payload": None,
            "expected": {
                "status": 200,
                "body_contains": ["orfs"],
                "each_orf_has": ["start", "end", "frame", "protein"],
            },
            "tags": ["analysis", "orfs"],
        })
        corpus["test_cases"].append({
            "id": f"analysis_composition_{gs['gene_symbol']}",
            "endpoint": "POST /analysis/composition",
            "query": {"sequence": gs["seq"]},
            "payload": None,
            "expected": {
                "status": 200,
                "body_contains": ["gc_content", "symbol_counts"],
                "gc_range": [0.0, 100.0],
            },
            "tags": ["analysis", "composition"],
        })

    # ── Protein properties ────────────────────────────────────────────────────
    corpus["test_cases"].append({
        "id": "analysis_protein_props_hbb",
        "endpoint": "POST /analysis/protein_properties",
        "query": {"sequence": MSA_SEQS[1]["seq"]},
        "payload": None,
        "expected": {
            "status": 200,
            "body_contains": ["molecular_weight_da", "isoelectric_point", "instability_index"],
            "mw_range": [5000, 200000],
            "pi_range": [3.0, 12.0],
        },
        "tags": ["analysis", "protein"],
    })

    # ── Database search ───────────────────────────────────────────────────────
    for symbol, _, desc in GENE_QUERIES[:3]:
        corpus["test_cases"].append({
            "id": f"search_{symbol}",
            "endpoint": "POST /databases/search",
            "payload": {
                "query": symbol,
                "databases": ["ncbi_gene", "ensembl"],
                "species": "human",
                "max_results": 5,
            },
            "expected": {
                "status": 200,
                "body_contains": ["results", "total"],
                "note": "If NCBI errors, results should NOT contain error text as a row",
                "each_result_has": ["id", "name", "database"],
            },
            "tags": ["databases", "search"],
        })

    print(f"\n  Corpus ready: {len(corpus['test_cases'])} test cases")
    return corpus


def main():
    corpus = build_corpus()
    out = CORPUS_DIR / "corpus.json"
    with open(out, "w") as f:
        json.dump(corpus, f, indent=2)
    print(f"\n✓ Corpus written to {out}")

    # Update agent_comms
    obj = _state.load()
    obj["workflow_state"]["corpus_ready"] = True
    obj["workflow_state"]["current_phase"] = 2
    obj["workflow_state"]["phases_completed"] = ["discovery", "analysis"]
    obj["messages"].append({
        "id": "msg_1",
        "from": "scout",
        "to": "runner",
        "type": "handoff",
        "subject": "Corpus ready — proceed to execution",
        "status": "pending",
        "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "body": (
            f"Discovery and analysis complete. {len(corpus['test_cases'])} test cases "
            f"written to .claude/qa/corpus/corpus.json. "
            f"Start backend with `python run_backend.py` then run runner.py."
        ),
    })
    _state.save(obj)
    print("✓ pipeline_state.json updated")


if __name__ == "__main__":
    main()
