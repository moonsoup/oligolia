"""gnomAD GraphQL API client — https://gnomad.broadinstitute.org/api"""

import httpx
from typing import Any

GRAPHQL_URL = "https://gnomad.broadinstitute.org/api"


class GnomADClient:
    def __init__(self, timeout: int = 30) -> None:
        self._timeout = timeout

    def _query(self, query: str, variables: dict | None = None) -> dict[str, Any]:
        payload: dict[str, Any] = {"query": query}
        if variables:
            payload["variables"] = variables
        with httpx.Client(timeout=self._timeout) as client:
            r = client.post(GRAPHQL_URL, json=payload,
                            headers={"Content-Type": "application/json"})
            r.raise_for_status()
            data = r.json()
            if "errors" in data:
                raise ValueError(f"gnomAD GraphQL errors: {data['errors']}")
            return data.get("data", {})

    def gene_by_symbol(self, symbol: str, dataset: str = "gnomad_r4") -> dict:
        query = """
        query GeneBySymbol($symbol: String!, $dataset: DatasetId!) {
          gene(gene_symbol: $symbol, reference_genome: GRCh38) {
            gene_id
            gene_version
            symbol
            hgnc_id
            name
            chrom
            start
            stop
            strand
            gnomad_constraint {
              exp_lof
              obs_lof
              lof_z
              pLI
            }
          }
        }
        """
        return self._query(query, {"symbol": symbol, "dataset": dataset})

    def variant(self, variant_id: str, dataset: str = "gnomad_r4") -> dict:
        """Fetch allele frequency and quality for a variant (e.g. 1-55516888-G-GA)."""
        query = """
        query Variant($variantId: String!, $dataset: DatasetId!) {
          variant(variant_id: $variantId, dataset: $dataset) {
            variant_id
            chrom
            pos
            ref
            alt
            genome {
              ac
              an
              af
              homozygote_count
            }
            exome {
              ac
              an
              af
              homozygote_count
            }
            rsids
            in_silico_predictors { id value }
          }
        }
        """
        return self._query(query, {"variantId": variant_id, "dataset": dataset})

    def variants_in_gene(self, gene_id: str, dataset: str = "gnomad_r4") -> dict:
        query = """
        query VariantsInGene($geneId: String!, $dataset: DatasetId!) {
          gene(gene_id: $geneId, reference_genome: GRCh38) {
            variants(dataset: $dataset) {
              variant_id
              pos
              ref
              alt
              genome { af }
              exome { af }
              consequence
            }
          }
        }
        """
        return self._query(query, {"geneId": gene_id, "dataset": dataset})

    def get_af(self, variant_id: str, dataset: str = "gnomad_r4") -> float | None:
        """Get global allele frequency for a variant, or None if not in gnomAD."""
        try:
            data = self.variant(variant_id, dataset)
            v = data.get("variant", {})
            genome_af = (v.get("genome") or {}).get("af")
            exome_af = (v.get("exome") or {}).get("af")
            if genome_af is not None:
                return genome_af
            return exome_af
        except Exception:
            return None
