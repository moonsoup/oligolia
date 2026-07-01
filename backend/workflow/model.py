"""Serializable workflow model — steps, status, and the workflow container."""

from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, Field


class StepType(str, Enum):
    """The composable operations a workflow step can run."""

    DB_SEARCH = "db_search"
    CRISPR_DESIGN = "crispr_design"
    OFF_TARGET = "off_target"
    PRIMER_DESIGN = "primer_design"
    CODON_OPTIMIZE = "codon_optimize"
    RESTRICTION_DIGEST = "restriction_digest"
    MSA = "msa"
    EXPORT = "export"
    ORDER = "order"


class StepStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETE = "complete"
    FAILED = "failed"
    SKIPPED = "skipped"


class WorkflowStep(BaseModel):
    """One operation in a workflow.

    ``params`` holds the step's own inputs; results and errors are populated at
    run time. ``result`` is stored as JSON-native data so the whole workflow
    round-trips cleanly through ``.ogo`` serialization.
    """

    id: str
    type: StepType
    name: str = ""
    params: dict[str, Any] = Field(default_factory=dict)
    status: StepStatus = StepStatus.PENDING
    error: str | None = None
    result: Any | None = None

    def reset(self) -> None:
        """Return the step to a pending, resultless state (for re-runs)."""
        self.status = StepStatus.PENDING
        self.error = None
        self.result = None


class Workflow(BaseModel):
    """A named, ordered pipeline of steps."""

    name: str
    description: str = ""
    format_version: str = "1.0"
    steps: list[WorkflowStep] = Field(default_factory=list)

    def step(self, step_id: str) -> WorkflowStep | None:
        return next((s for s in self.steps if s.id == step_id), None)

    def reset(self) -> None:
        for s in self.steps:
            s.reset()

    @property
    def is_complete(self) -> bool:
        return bool(self.steps) and all(
            s.status == StepStatus.COMPLETE for s in self.steps
        )
