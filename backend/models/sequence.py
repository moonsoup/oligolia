from enum import Enum
from typing import Any
from pydantic import BaseModel, Field


class MoleculeType(str, Enum):
    DNA = "DNA"
    RNA = "RNA"
    PROTEIN = "PROTEIN"
    UNKNOWN = "UNKNOWN"


class Strand(str, Enum):
    PLUS = "+"
    MINUS = "-"
    BOTH = "."


class Annotation(BaseModel):
    feature_type: str
    start: int
    end: int
    strand: Strand = Strand.PLUS
    qualifiers: dict[str, Any] = Field(default_factory=dict)


class Sequence(BaseModel):
    id: str
    name: str = ""
    description: str = ""
    seq: str
    molecule_type: MoleculeType = MoleculeType.DNA
    annotations: list[Annotation] = Field(default_factory=list)
    source_db: str = ""
    accession: str = ""
    length: int = 0

    def model_post_init(self, __context: Any) -> None:
        if not self.length:
            self.length = len(self.seq)


class SequenceEditRequest(BaseModel):
    operation: str  # insert | delete | replace | reverse_complement | translate | complement
    position: int | None = None
    end_position: int | None = None
    insert_seq: str | None = None
    replacement: str | None = None


class SequenceEditResult(BaseModel):
    original_id: str
    operation: str
    result_seq: str
    diff_start: int | None = None
    diff_end: int | None = None
    message: str = ""
