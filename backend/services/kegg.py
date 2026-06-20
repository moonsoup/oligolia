"""KEGG REST API client — https://www.kegg.jp/kegg/rest/keggapi.html"""

import httpx
from typing import Any

BASE = "https://rest.kegg.jp"


class KEGGClient:
    def __init__(self, timeout: int = 30) -> None:
        self._timeout = timeout

    def _get(self, path: str) -> str:
        with httpx.Client(timeout=self._timeout) as client:
            r = client.get(f"{BASE}{path}")
            r.raise_for_status()
            return r.text

    def _get_json(self, path: str) -> Any:
        with httpx.Client(timeout=self._timeout) as client:
            r = client.get(f"{BASE}{path}", headers={"Accept": "application/json"})
            r.raise_for_status()
            return r.json()

    def search(self, database: str, query: str) -> list[tuple[str, str]]:
        """Search KEGG database. Returns list of (kegg_id, description) tuples."""
        text = self._get(f"/find/{database}/{query}")
        results = []
        for line in text.strip().splitlines():
            if "\t" in line:
                parts = line.split("\t", 1)
                results.append((parts[0], parts[1]))
        return results

    def get_entry(self, entry_id: str) -> str:
        """Fetch a KEGG entry (gene, pathway, compound, etc.)."""
        return self._get(f"/get/{entry_id}")

    def get_pathway(self, pathway_id: str) -> str:
        """Fetch pathway record (e.g. hsa05210 = colorectal cancer)."""
        return self._get(f"/get/{pathway_id}")

    def get_pathway_image(self, pathway_id: str) -> bytes:
        """Fetch pathway PNG image."""
        with httpx.Client(timeout=self._timeout) as client:
            r = client.get(f"{BASE}/get/{pathway_id}/image")
            r.raise_for_status()
            return r.content

    def list_pathways(self, organism: str = "hsa") -> list[tuple[str, str]]:
        """List all pathways for an organism (hsa = human)."""
        return self.search("pathway", organism)

    def genes_in_pathway(self, pathway_id: str) -> list[str]:
        """Return KEGG gene IDs for all genes in a pathway."""
        entry = self.get_entry(pathway_id)
        genes = []
        in_gene_section = False
        for line in entry.splitlines():
            if line.startswith("GENE"):
                in_gene_section = True
            elif in_gene_section:
                if line.startswith(" "):
                    parts = line.strip().split()
                    if parts:
                        genes.append(parts[0])
                else:
                    break
        return genes

    def convert(self, from_db: str, to_db: str, ids: list[str]) -> dict[str, str]:
        """Convert between identifier systems (e.g. ncbi-geneid → kegg genes)."""
        id_str = "+".join(ids)
        text = self._get(f"/conv/{to_db}/{from_db}:{id_str}")
        mapping: dict[str, str] = {}
        for line in text.strip().splitlines():
            if "\t" in line:
                src, dst = line.split("\t", 1)
                mapping[src] = dst
        return mapping

    def find_genes_by_name(self, gene_name: str, organism: str = "hsa") -> list[tuple[str, str]]:
        return self.search(f"genes/{organism}", gene_name)
