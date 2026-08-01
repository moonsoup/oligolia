from enum import Enum
from pydantic import BaseModel, Field


class StructureSource(str, Enum):
    EXPERIMENTAL_PDB = "experimental_pdb"
    PREDICTED_ALPHAFOLD_DB = "predicted_alphafold_db"
    PREDICTED_ESMFOLD = "predicted_esmfold"


class StructureRequest(BaseModel):
    sequence: str
    gene_symbol: str | None = None
    uniprot_id: str | None = None
    skip_experimental_lookup: bool = False


class StructureResult(BaseModel):
    source: StructureSource
    pdb_id: str | None = None
    pdb_text: str
    sequence_length: int
    confidence_note: str = ""
    warnings: list[str] = Field(default_factory=list)


class InteractionPoint(BaseModel):
    residue_index: int
    residue_name: str
    chain: str = "A"
    classification: str  # acidic | basic | polar | hydrophobic
    relative_sasa: float
    is_putative_interaction_point: bool


class InteractionPointsRequest(BaseModel):
    pdb_text: str


class InteractionPointsResult(BaseModel):
    points: list[InteractionPoint]
    method: str = (
        "Charge/polarity classification (static side-chain lookup table) + "
        "Shrake-Rupley relative solvent-accessible surface area, thresholded. "
        "Heuristic surface/charge annotation — NOT a validated binding-site or "
        "docking prediction."
    )
