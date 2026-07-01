from pydantic import BaseModel, Field
from enum import Enum


class CasType(str, Enum):
    CAS9 = "SpCas9"
    CAS9_HF = "SpCas9-HF1"
    CAS12A = "AsCas12a"
    CAS13 = "LwaCas13a"


class GuideRNA(BaseModel):
    sequence: str  # 20 nt protospacer
    pam: str
    position: int
    strand: str
    gc_content: float
    on_target_score: float | None = None
    off_target_count: int | None = None
    # Mismatch-bucketed off-target hits, e.g. {"0": 0, "1": 3, "2": 47}.
    off_target_summary: dict[str, int] | None = None
    # MIT-style specificity score (0–100); higher = fewer/weaker off-targets.
    specificity_score: float | None = None
    efficiency_score: float | None = None


class CRISPRDesignRequest(BaseModel):
    target_sequence: str
    cas_type: CasType = CasType.CAS9
    guide_length: int = Field(default=20, ge=17, le=24)
    max_guides: int = Field(default=10, ge=1, le=50)
    check_off_targets: bool = False
    # Extra reference loci to scan for off-targets. The target sequence is
    # always scanned; these are additional (e.g. paralogs or genomic context).
    reference_sequences: list[str] = Field(default_factory=list)
    max_mismatches: int = Field(default=3, ge=1, le=4)


class CRISPRDesignResponse(BaseModel):
    target_length: int
    cas_type: CasType
    guides: list[GuideRNA]
    total_candidates: int
