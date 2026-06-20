"""STRING protein interaction API — https://string-db.org/cgi/help.pl?subpage=api"""

import httpx
from typing import Any

BASE = "https://string-db.org/api"


class STRINGClient:
    def __init__(self, version: str = "12.0", timeout: int = 30) -> None:
        self.version = version
        self._timeout = timeout

    def _get(self, endpoint: str, params: dict) -> Any:
        url = f"{BASE}/json/{endpoint}"
        params["version"] = self.version
        params.setdefault("caller_identity", "oligolia_gene_editor")
        with httpx.Client(timeout=self._timeout) as client:
            r = client.get(url, params=params)
            r.raise_for_status()
            return r.json()

    def get_string_ids(self, identifiers: list[str], species: int = 9606) -> list[dict]:
        """Map gene names/identifiers to STRING IDs (9606 = human)."""
        return self._get("get_string_ids", {
            "identifiers": "\r".join(identifiers),
            "species": species,
        })

    def network(self, identifiers: list[str], species: int = 9606,
                required_score: int = 400, network_type: str = "functional") -> list[dict]:
        """Get protein-protein interaction network edges."""
        return self._get("network", {
            "identifiers": "\r".join(identifiers),
            "species": species,
            "required_score": required_score,
            "network_type": network_type,
        })

    def interaction_partners(self, identifier: str, species: int = 9606,
                             limit: int = 10, required_score: int = 700) -> list[dict]:
        """Get top interaction partners for a single protein."""
        return self._get("interaction_partners", {
            "identifier": identifier,
            "species": species,
            "limit": limit,
            "required_score": required_score,
        })

    def enrichment(self, identifiers: list[str], species: int = 9606) -> list[dict]:
        """GO/KEGG/Reactome enrichment for a set of proteins."""
        return self._get("enrichment", {
            "identifiers": "\r".join(identifiers),
            "species": species,
        })

    def get_image_url(self, identifiers: list[str], species: int = 9606) -> str:
        """Return URL to network image PNG."""
        params = {
            "identifiers": "%0d".join(identifiers),
            "species": species,
            "version": self.version,
        }
        query = "&".join(f"{k}={v}" for k, v in params.items())
        return f"{BASE}/image/network?{query}"
