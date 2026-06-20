"""NCBI Entrez E-utilities + Datasets v2 REST API client."""

import time
import httpx
from typing import Any

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
# Rate limit: 3 req/s without API key, 10 req/s with key
_last_request: float = 0.0


def _throttle() -> None:
    global _last_request
    elapsed = time.monotonic() - _last_request
    if elapsed < 0.34:
        time.sleep(0.34 - elapsed)
    _last_request = time.monotonic()


class NCBIClient:
    def __init__(self, api_key: str = "", email: str = "") -> None:
        self.api_key = api_key
        self.email = email
        self._params: dict[str, str] = {}
        if api_key:
            self._params["api_key"] = api_key
        if email:
            self._params["email"] = email

    def _get(self, url: str, params: dict) -> httpx.Response:
        _throttle()
        merged = {**self._params, **params}
        with httpx.Client(timeout=30) as client:
            r = client.get(url, params=merged)
            r.raise_for_status()
            return r

    # ── ESearch ──────────────────────────────────────────────────────────────

    def esearch(self, db: str, term: str, retmax: int = 20) -> dict[str, Any]:
        """Return list of UIDs matching term in db."""
        r = self._get(
            f"{ENTREZ_BASE}/esearch.fcgi",
            {"db": db, "term": term, "retmax": retmax, "retmode": "json"},
        )
        return r.json()

    # ── EFetch ───────────────────────────────────────────────────────────────

    def efetch(self, db: str, ids: list[str], rettype: str = "gb", retmode: str = "text") -> str:
        """Fetch full records for given UIDs."""
        r = self._get(
            f"{ENTREZ_BASE}/efetch.fcgi",
            {"db": db, "id": ",".join(ids), "rettype": rettype, "retmode": retmode},
        )
        return r.text

    def efetch_json(self, db: str, ids: list[str]) -> dict[str, Any]:
        """Fetch records in JSON summary format."""
        r = self._get(
            f"{ENTREZ_BASE}/esummary.fcgi",
            {"db": db, "id": ",".join(ids), "retmode": "json"},
        )
        return r.json()

    # ── ELink ────────────────────────────────────────────────────────────────

    def elink(self, dbfrom: str, db: str, ids: list[str]) -> dict[str, Any]:
        r = self._get(
            f"{ENTREZ_BASE}/elink.fcgi",
            {"dbfrom": dbfrom, "db": db, "id": ",".join(ids), "retmode": "json"},
        )
        return r.json()

    # ── Convenience ──────────────────────────────────────────────────────────

    def search_genes(self, query: str, organism: str = "Homo sapiens", max_results: int = 20) -> list[dict]:
        term = f"{query}[Gene Name] AND {organism}[Organism]"
        result = self.esearch("gene", term, retmax=max_results)
        ids = result.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        summary = self.efetch_json("gene", ids)
        records = summary.get("result", {})
        uids = records.get("uids", [])
        return [records[uid] for uid in uids if uid in records]

    def fetch_nucleotide(self, accession: str) -> str:
        """Fetch GenBank flat file for a nucleotide accession."""
        return self.efetch("nucleotide", [accession], rettype="gb", retmode="text")

    def fetch_fasta(self, accession: str, db: str = "nucleotide") -> str:
        return self.efetch(db, [accession], rettype="fasta", retmode="text")

    def blast_search(self, sequence: str, program: str = "blastn", database: str = "nt", max_hits: int = 10) -> dict:
        """Submit a BLAST search via NCBI BLAST REST API and poll for results."""
        put_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
        with httpx.Client(timeout=60) as client:
            r = client.post(put_url, data={
                "CMD": "Put",
                "PROGRAM": program,
                "DATABASE": database,
                "QUERY": sequence,
                "FORMAT_TYPE": "JSON2",
                "HITLIST_SIZE": max_hits,
            })
            r.raise_for_status()
            # Extract RID
            rid = None
            for line in r.text.splitlines():
                if line.startswith("    RID = "):
                    rid = line.split("=")[1].strip()
                    break
            if not rid:
                raise ValueError("BLAST: could not obtain RID from PUT response")

            # Poll for results
            for _ in range(30):
                time.sleep(5)
                status_r = client.get(put_url, params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2"})
                if "Status=WAITING" in status_r.text:
                    continue
                if "Status=READY" in status_r.text or status_r.headers.get("content-type", "").startswith("application/json"):
                    try:
                        return status_r.json()
                    except Exception:
                        # Try extracting JSON from mixed response
                        start = status_r.text.find("{")
                        if start >= 0:
                            import json
                            return json.loads(status_r.text[start:])
                if "Status=FAILED" in status_r.text:
                    raise RuntimeError("BLAST search failed")
            raise TimeoutError("BLAST search timed out after 150 seconds")

    def search_clinvar(self, gene: str, max_results: int = 50) -> list[dict]:
        term = f"{gene}[Gene Name]"
        result = self.esearch("clinvar", term, retmax=max_results)
        ids = result.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        summary = self.efetch_json("clinvar", ids)
        records = summary.get("result", {})
        uids = records.get("uids", [])
        return [records[uid] for uid in uids if uid in records]
