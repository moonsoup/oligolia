"""AlphaFold Protein Structure Database API client — https://alphafold.ebi.ac.uk/api-docs

Free, keyless lookup of precomputed structure predictions covering essentially the
whole UniProt proteome. Preferred over live-folding (esmfold.py) whenever a UniProt
accession is known: no length cap, no compute wait, confidence already computed.
"""

import httpx
from pydantic import BaseModel

BASE = "https://alphafold.ebi.ac.uk/api/prediction"


class AlphaFoldPrediction(BaseModel):
    uniprot_accession: str
    pdb_text: str
    model_created_date: str = ""
    global_metric_value: float | None = None  # overall pLDDT
    fraction_plddt_very_low: float | None = None
    fraction_plddt_low: float | None = None
    fraction_plddt_confident: float | None = None
    fraction_plddt_very_high: float | None = None


class AlphaFoldDBClient:
    def __init__(self, timeout: int = 30) -> None:
        self._timeout = timeout

    def get_prediction(self, uniprot_accession: str) -> AlphaFoldPrediction | None:
        """Fetch the precomputed AlphaFold DB entry for a UniProt accession.

        Returns None if AlphaFold DB has no model for this accession (a 404 —
        rare, but happens for very short or otherwise excluded entries).
        """
        with httpx.Client(timeout=self._timeout) as client:
            r = client.get(f"{BASE}/{uniprot_accession}")
            if r.status_code == 404:
                return None
            r.raise_for_status()
            entries = r.json()
            if not entries:
                return None
            entry = entries[0]

            pdb_url = entry.get("pdbUrl")
            if not pdb_url:
                return None
            pdb_r = client.get(pdb_url)
            pdb_r.raise_for_status()

            return AlphaFoldPrediction(
                uniprot_accession=uniprot_accession,
                pdb_text=pdb_r.text,
                model_created_date=entry.get("modelCreatedDate", ""),
                global_metric_value=entry.get("globalMetricValue"),
                fraction_plddt_very_low=entry.get("fractionPlddtVeryLow"),
                fraction_plddt_low=entry.get("fractionPlddtLow"),
                fraction_plddt_confident=entry.get("fractionPlddtConfident"),
                fraction_plddt_very_high=entry.get("fractionPlddtVeryHigh"),
            )
