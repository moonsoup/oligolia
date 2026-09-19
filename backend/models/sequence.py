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


class LocationPart(BaseModel):
    """Everything an INSDC location part carries besides its two integers (#94).

    Index-aligned with `Annotation.parts`. Flattening a Biopython location to
    `(int, int)` pairs threw all of this away before it was ever stored, so
    export could not recover it: a `<108..1007` CDS came back as `108..1007`,
    which asserts a start codon the original explicitly declined to claim.

    `start_class` / `end_class` are the position *class* of each boundary:
    `"before"` is INSDC's `<` ("the feature starts before the first sequenced
    base"), `"after"` its `>`, `"exact"` a plain coordinate. They cannot be
    inferred from the numbers — Biopython's `BeforePosition` subclasses `int`
    and `BeforePosition(5) == ExactPosition(5)`, so no coordinate comparison can
    see the difference.

    `ref` is the remote accession of a part like `J00194.1:100..202`, whose
    bases belong to another entry entirely. Dropping it does not merely lose
    information, it reassigns the feature to this record's own bases 100..202.
    """

    start_class: str = "exact"
    end_class: str = "exact"
    ref: str | None = None


class Annotation(BaseModel):
    feature_type: str
    #: Outer bounds, 0-based half-open. For a spliced or origin-spanning feature
    #: these are the min/max across `parts`, which is what every existing
    #: consumer (plasmid map, feature table, GUI) already reads — so their
    #: meaning is unchanged.
    start: int
    end: int
    strand: Strand = Strand.PLUS
    #: Every interval the feature actually occupies, in file order. A simple
    #: feature has exactly one part equal to (start, end); a GenBank
    #: `join(5..10,30..40)` has two.
    #:
    #: Added for #57: keeping only start/end collapsed every join to its outer
    #: bounds, so an intron was swallowed and `join(91..100,1..20)` on a 100 bp
    #: circular record loaded as the whole plasmid — which the plasmid map then
    #: drew as a full circle.
    parts: list[tuple[int, int]] = Field(default_factory=list)
    #: Per-part location detail, index-aligned with `parts` (#94). Empty means
    #: "every boundary exact, no remote reference", which is what a reader with
    #: no notion of partial boundaries (GFF, FASTA, the cloning planner) gives —
    #: so it stays the default and every existing producer keeps working.
    part_details: list[LocationPart] = Field(default_factory=list)
    #: The INSDC compound operator, `"join"` or `"order"` (#94). Only meaningful
    #: with more than one part. `join` asserts the parts form one contiguous
    #: sequence; `order` asserts only that they occur in this order and says
    #: nothing about joining them, so writing one as the other changes the claim.
    location_operator: str = "join"
    qualifiers: dict[str, Any] = Field(default_factory=dict)

    def details_per_part(self) -> list["LocationPart"]:
        """`part_details`, padded to exactly one entry per part (#94).

        Callers get a same-length list whether or not the source carried any
        detail, so nothing has to special-case the empty default.
        """
        count = len(self.parts) or 1
        details = list(self.part_details[:count])
        details.extend(LocationPart() for _ in range(count - len(details)))
        return details


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
    is_circular: bool = False

    def model_post_init(self, __context: Any) -> None:
        if not self.length:
            self.length = len(self.seq)


class ReferenceFeature(BaseModel):
    """A known vector part in the bundled reference library (issue #42).

    ``translation`` is set for coding parts (the peptide the DNA should encode)
    and None for non-coding parts (promoters, operators, terminators, primer
    binding sites). The homology-scan feature (#43) consumes these records.
    """

    name: str
    feature_type: str
    sequence: str
    molecule_type: MoleculeType = MoleculeType.DNA
    translation: str | None = None
    source: str = ""


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
