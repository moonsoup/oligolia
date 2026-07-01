"""Multi-step workflow builder — chain Oligolia operations into a tracked pipeline.

This is the backend foundation for issue #10: a serializable workflow model,
a pure execution engine that runs steps sequentially (output of step N feeds
step N+1 via a shared context), and ``.ogo`` save/load. The GUI Workflow tab
and any live ordering-API submission are separate follow-ups; the ``order``
step here reuses the offline vendor file export (no order is submitted).
"""

from .model import StepStatus, StepType, Workflow, WorkflowStep
from .engine import run_workflow, run_step, register_handler, HANDLERS
from .ogo import dumps, loads, save_ogo, load_ogo

__all__ = [
    "StepStatus", "StepType", "Workflow", "WorkflowStep",
    "run_workflow", "run_step", "register_handler", "HANDLERS",
    "dumps", "loads", "save_ogo", "load_ogo",
]
