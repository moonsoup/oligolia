from pydantic import BaseModel, Field
from enum import Enum


class VariantType(str, Enum):
    SNP = "SNP"
    INDEL = "INDEL"
    DEL = "DEL"
    INS = "INS"
    CNV = "CNV"
    SV = "SV"
    UNKNOWN = "UNKNOWN"


class ClinicalSignificance(str, Enum):
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely pathogenic"
    VUS = "Uncertain significance"
    LIKELY_BENIGN = "Likely benign"
    BENIGN = "Benign"
    CONFLICTING = "Conflicting interpretations"
    NOT_PROVIDED = "Not provided"


class Variant(BaseModel):
    chrom: str
    pos: int
    ref: str
    alt: list[str]
    id: str = "."
    qual: float | None = None
    filter: list[str] = Field(default_factory=list)
    info: dict = Field(default_factory=dict)
    variant_type: VariantType = VariantType.UNKNOWN
    gene: str = ""
    clinical_significance: ClinicalSignificance | None = None
    allele_frequency: float | None = None
    gnomad_af: float | None = None


class VariantAnnotationRequest(BaseModel):
    variants: list[Variant]
    annotate_clinvar: bool = True
    annotate_gnomad: bool = True


class VariantAnnotationResponse(BaseModel):
    variants: list[Variant]
    total: int
    annotated: int
