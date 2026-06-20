"""Ensembl REST API client — https://rest.ensembl.org"""

import httpx
from typing import Any

BASE = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}


class EnsemblClient:
    def __init__(self, timeout: int = 30) -> None:
        self._timeout = timeout

    def _get(self, endpoint: str, params: dict | None = None) -> Any:
        url = f"{BASE}{endpoint}"
        with httpx.Client(timeout=self._timeout, headers=HEADERS) as client:
            r = client.get(url, params=params or {})
            r.raise_for_status()
            return r.json()

    def _post(self, endpoint: str, body: dict) -> Any:
        url = f"{BASE}{endpoint}"
        with httpx.Client(timeout=self._timeout, headers=HEADERS) as client:
            r = client.post(url, json=body)
            r.raise_for_status()
            return r.json()

    # ── Gene lookup ──────────────────────────────────────────────────────────

    def lookup_symbol(self, species: str, symbol: str) -> dict:
        """Look up a gene by symbol (e.g. BRCA2, TP53)."""
        return self._get(f"/lookup/symbol/{species}/{symbol}", {"expand": 1})

    def lookup_id(self, ensembl_id: str) -> dict:
        """Look up any feature by Ensembl stable ID."""
        return self._get(f"/lookup/id/{ensembl_id}", {"expand": 1})

    def search(self, species: str, query: str) -> list[dict]:
        """Free-text search for genes within a species."""
        result = self._get(f"/lookup/symbol/{species}/{query}")
        return [result] if result else []

    # ── Sequence ─────────────────────────────────────────────────────────────

    def sequence_id(self, ensembl_id: str, seq_type: str = "genomic") -> dict:
        """Fetch sequence for an Ensembl ID (genomic, cds, cdna, protein)."""
        return self._get(f"/sequence/id/{ensembl_id}", {"type": seq_type})

    def sequence_region(self, species: str, region: str) -> dict:
        """Fetch sequence for a genomic region (e.g. 'X:1000000..1000100')."""
        return self._get(f"/sequence/region/{species}/{region}")

    # ── Variants ─────────────────────────────────────────────────────────────

    def variation(self, species: str, rsid: str) -> dict:
        """Get variant details by rsID."""
        return self._get(f"/variation/{species}/{rsid}")

    def vep_hgvs(self, species: str, hgvs: str) -> list[dict]:
        """Variant Effect Predictor for an HGVS notation."""
        return self._get(f"/vep/{species}/hgvs/{hgvs}")

    def vep_region(self, species: str, region: str, allele: str) -> list[dict]:
        """VEP for a region + allele string."""
        return self._get(f"/vep/{species}/region/{region}/{allele}")

    # ── Homology ─────────────────────────────────────────────────────────────

    def homology_id(self, species: str, ensembl_id: str, target_species: str | None = None) -> dict:
        params: dict = {"type": "all"}
        if target_species:
            params["target_species"] = target_species
        return self._get(f"/homology/id/{species}/{ensembl_id}", params)

    # ── Phenotype ────────────────────────────────────────────────────────────

    def phenotype_gene(self, species: str, gene: str) -> list[dict]:
        return self._get(f"/phenotype/gene/{species}/{gene}")

    # ── Features on region ───────────────────────────────────────────────────

    def overlap_region(self, species: str, region: str, feature: str = "gene") -> list[dict]:
        return self._get(f"/overlap/region/{species}/{region}", {"feature": feature})

    # ── Cross-references ─────────────────────────────────────────────────────

    def xrefs_id(self, ensembl_id: str) -> list[dict]:
        return self._get(f"/xrefs/id/{ensembl_id}")

    def xrefs_symbol(self, species: str, symbol: str) -> list[dict]:
        return self._get(f"/xrefs/symbol/{species}/{symbol}")
