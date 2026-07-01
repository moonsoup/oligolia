"""``.ogo`` workflow serialization — a plain-JSON save/load format.

An ``.ogo`` file is the JSON dump of a :class:`Workflow`, including each step's
last-run status and result, so a saved workflow can be shared, re-opened, and
re-run on a new gene.
"""

from __future__ import annotations

from pathlib import Path

from .model import Workflow


def dumps(wf: Workflow) -> str:
    """Serialize a workflow to a pretty-printed JSON string."""
    return wf.model_dump_json(indent=2)


def loads(text: str) -> Workflow:
    """Parse a workflow from an ``.ogo`` JSON string."""
    return Workflow.model_validate_json(text)


def save_ogo(wf: Workflow, path: str | Path) -> None:
    Path(path).write_text(dumps(wf), encoding="utf-8")


def load_ogo(path: str | Path) -> Workflow:
    return loads(Path(path).read_text(encoding="utf-8"))
