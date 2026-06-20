"""UniProt REST API client — https://www.uniprot.org/help/uniprot_rest_tutorial"""

import httpx
from typing import Any

BASE = "https://rest.uniprot.org"


class UniProtClient:
    def __init__(self, timeout: int = 30) -> None:
        self._timeout = timeout

    def _get(self, path: str, params: dict | None = None) -> Any:
        with httpx.Client(timeout=self._timeout) as client:
            r = client.get(f"{BASE}{path}", params=params or {})
            r.raise_for_status()
            ct = r.headers.get("content-type", "")
            return r.json() if "json" in ct else r.text

    def search(self, query: str, fields: str = "accession,id,protein_name,gene_names,organism_name,length",
               format: str = "json", size: int = 25) -> dict:
        """Search UniProtKB."""
        return self._get("/uniprotkb/search", {
            "query": query, "fields": fields, "format": format, "size": size
        })

    def get_entry(self, accession: str, format: str = "json") -> Any:
        """Fetch a single UniProt entry by accession (e.g. P04637 = TP53)."""
        return self._get(f"/uniprotkb/{accession}", {"format": format})

    def get_fasta(self, accession: str) -> str:
        return self._get(f"/uniprotkb/{accession}.fasta")

    def id_mapping(self, from_db: str, to_db: str, ids: list[str]) -> dict:
        """Map identifiers between databases (e.g. UniProtKB_AC → Ensembl)."""
        with httpx.Client(timeout=60) as client:
            r = client.post(
                f"{BASE}/idmapping/run",
                data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
            )
            r.raise_for_status()
            job_id = r.json()["jobId"]

            import time
            for _ in range(20):
                time.sleep(2)
                status = client.get(f"{BASE}/idmapping/status/{job_id}")
                status.raise_for_status()
                data = status.json()
                if data.get("jobStatus") == "FINISHED" or "results" in data:
                    results_r = client.get(f"{BASE}/idmapping/results/{job_id}")
                    results_r.raise_for_status()
                    return results_r.json()
            raise TimeoutError("UniProt ID mapping timed out")

    def search_by_gene(self, gene_name: str, organism: str = "homo sapiens", size: int = 10) -> dict:
        query = f"gene:{gene_name} AND organism_name:{organism} AND reviewed:true"
        return self.search(query, size=size)
