"""Reactome Content Service API — https://reactome.org/dev/content-service"""

import httpx
from typing import Any

BASE = "https://reactome.org/ContentService"


class ReactomeClient:
    def __init__(self, timeout: int = 30) -> None:
        self._timeout = timeout

    def _get(self, path: str, params: dict | None = None) -> Any:
        with httpx.Client(timeout=self._timeout) as client:
            r = client.get(f"{BASE}{path}", params=params or {},
                           headers={"Accept": "application/json"})
            r.raise_for_status()
            return r.json()

    def search_query(self, query: str, species: str = "Homo sapiens", page: int = 1, page_size: int = 20) -> dict:
        """Full-text search across Reactome."""
        return self._get("/search/query", {
            "query": query, "species": species, "page": page, "pageSize": page_size
        })

    def pathways_for_gene(self, gene_id: str, species: str = "Homo sapiens") -> list[dict]:
        """Get pathways that contain a gene (by Ensembl or UniProt ID)."""
        return self._get(f"/data/pathways/low/entity/{gene_id}", {"species": species})

    def pathway_hierarchy(self, species: str = "Homo sapiens") -> list[dict]:
        return self._get(f"/data/eventsHierarchy/{species}")

    def get_pathway(self, pathway_id: str) -> dict:
        return self._get(f"/data/query/{pathway_id}")

    def pathway_participants(self, pathway_id: str) -> list[dict]:
        return self._get(f"/data/participants/{pathway_id}")

    def disease_pathways(self) -> list[dict]:
        """List all disease-related pathways."""
        return self._get("/data/diseases/doid")

    def mapping_analysis(self, identifiers: list[str], species: str = "Homo sapiens") -> dict:
        """Perform pathway enrichment analysis on a list of identifiers."""
        with httpx.Client(timeout=60) as client:
            r = client.post(
                f"{BASE}/analysis/identifiers",
                content="\n".join(identifiers),
                params={"species": species},
                headers={"Content-Type": "text/plain", "Accept": "application/json"},
            )
            r.raise_for_status()
            return r.json()
