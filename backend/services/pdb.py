"""RCSB PDB REST API client — https://data.rcsb.org/redoc/"""

import httpx
from typing import Any

BASE = "https://data.rcsb.org"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


class PDBClient:
    def __init__(self, timeout: int = 30) -> None:
        self._timeout = timeout

    def _get(self, path: str) -> Any:
        with httpx.Client(timeout=self._timeout) as client:
            r = client.get(f"{BASE}{path}")
            r.raise_for_status()
            return r.json()

    def get_entry(self, pdb_id: str) -> dict:
        """Fetch full metadata for a PDB entry."""
        return self._get(f"/rest/v1/core/entry/{pdb_id.upper()}")

    def get_polymer_entity(self, pdb_id: str, entity_id: int = 1) -> dict:
        return self._get(f"/rest/v1/core/polymer_entity/{pdb_id.upper()}/{entity_id}")

    def search_by_gene(self, gene_name: str, max_results: int = 10) -> list[str]:
        """Text search for PDB entries by gene name. Returns list of PDB IDs."""
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                    "operator": "exact_match",
                    "value": gene_name.upper(),
                },
            },
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": max_results}},
        }
        with httpx.Client(timeout=self._timeout) as client:
            r = client.post(SEARCH_URL, json=query)
            r.raise_for_status()
            data = r.json()
            return [hit["identifier"] for hit in data.get("result_set", [])]

    def search_by_uniprot(self, uniprot_id: str, max_results: int = 10) -> list[str]:
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": uniprot_id,
                },
            },
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": max_results}},
        }
        with httpx.Client(timeout=self._timeout) as client:
            r = client.post(SEARCH_URL, json=query)
            r.raise_for_status()
            data = r.json()
            return [hit["identifier"] for hit in data.get("result_set", [])]

    def download_pdb(self, pdb_id: str) -> str:
        """Download PDB flat file content."""
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        with httpx.Client(timeout=60) as client:
            r = client.get(url)
            r.raise_for_status()
            return r.text

    def download_mmcif(self, pdb_id: str) -> str:
        """Download mmCIF format."""
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        with httpx.Client(timeout=60) as client:
            r = client.get(url)
            r.raise_for_status()
            return r.text
