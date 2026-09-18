"""The QA pipeline's own inter-phase state — one definition, one file.

Each of the four phases used to build its own path:

    comms_path = BASE.parent / "agent_comms.json"   # -> .claude/agent_comms.json

`BASE.parent` is `.claude/`, not the repo root, so that file never existed and
every `if comms_path.exists():` block was dead — the phases never handed off and
the "updated" line never printed, which is why nothing looked wrong (#80).

Repointing it at the repo's `agent_comms.json` would have been worse: that file is
the Claude<->Codex review protocol (`_protocol` v1.3 + typed messages), and the
pipeline writes `workflow_state`, so the two would collide on one filename with
incompatible schemas. They are separate concerns and now have separate files.

Ziggurat flagged the same thing structurally: one path, four namers.
"""

from __future__ import annotations

import json
from pathlib import Path

#: Beside the scripts that own it, not in .claude/ and not in the repo root.
STATE_PATH = Path(__file__).parent / "pipeline_state.json"

PHASES = ["discover", "corpus", "execute", "compare", "locate", "report"]


def _empty() -> dict:
    return {
        "_about": (
            "Inter-phase state for the QA pipeline (.claude/qa). Not the "
            "Claude<->Codex review protocol — that is agent_comms.json in the "
            "repo root, with a different schema (#80)."
        ),
        "workflow_state": {"current_phase": 0, "phases_completed": []},
        "messages": [],
    }


def load() -> dict:
    """The current state, creating it if this is the first run."""
    if not STATE_PATH.exists():
        return _empty()
    try:
        state = json.loads(STATE_PATH.read_text())
    except (OSError, ValueError):
        return _empty()
    # Tolerate a file written by an older version.
    state.setdefault("workflow_state", {"current_phase": 0, "phases_completed": []})
    state["workflow_state"].setdefault("current_phase", 0)
    state["workflow_state"].setdefault("phases_completed", [])
    state.setdefault("messages", [])
    return state


def save(state: dict) -> None:
    STATE_PATH.write_text(json.dumps(state, indent=2) + "\n")


def record(phase: int, completed: list[str], message: dict | None = None) -> dict:
    """Mark phases done and optionally append a handoff message."""
    state = load()
    state["workflow_state"]["current_phase"] = phase
    for name in completed:
        if name not in state["workflow_state"]["phases_completed"]:
            state["workflow_state"]["phases_completed"].append(name)
    if message is not None:
        state["messages"].append(message)
    save(state)
    return state
