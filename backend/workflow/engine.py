"""Sequential workflow execution engine.

Each step's handler reads inputs from ``step.params`` and a shared ``context``
dict, and writes its outputs back into the context so the next step can consume
them (output of step N feeds step N+1). Handlers reuse the existing operation
functions; results are stored JSON-native so a run round-trips through ``.ogo``.

Handlers are registered per :class:`StepType`. The offline computational steps
are wired here; ``db_search``, ``codon_optimize`` and ``msa`` are intentionally
left unregistered in this backend slice and report a clear, non-fatal error.
"""

from __future__ import annotations

from typing import Any, Callable

from fastapi import HTTPException

from .model import StepStatus, StepType, Workflow, WorkflowStep

Handler = Callable[[WorkflowStep, dict], Any]
HANDLERS: dict[StepType, Handler] = {}


def register_handler(step_type: StepType, handler: Handler) -> None:
    """Register (or override) the handler for a step type."""
    HANDLERS[step_type] = handler


def _require_sequence(step: WorkflowStep, ctx: dict) -> str:
    seq = (
        step.params.get("sequence")
        or step.params.get("target_sequence")
        or step.params.get("template")
        or ctx.get("sequence")
    )
    if not seq:
        raise ValueError(
            "No sequence available — set params.sequence or run a step that "
            "produces one (e.g. db_search) earlier in the workflow."
        )
    return seq


# ── Handlers ─────────────────────────────────────────────────────────────────

def _h_crispr_design(step: WorkflowStep, ctx: dict) -> dict:
    from ..models.crispr import CRISPRDesignRequest
    from ..routers.crispr import design_guides

    seq = _require_sequence(step, ctx)
    payload = dict(step.params)
    payload.pop("sequence", None)
    payload.setdefault("target_sequence", seq)
    resp = design_guides(CRISPRDesignRequest(**payload))
    ctx["sequence"] = seq
    ctx["guides"] = [g.model_dump() for g in resp.guides]
    return resp.model_dump()


def _h_off_target(step: WorkflowStep, ctx: dict) -> dict:
    from ..crispr_offtarget import scan_off_targets

    guides = step.params.get("guides") or ctx.get("guides")
    if not guides:
        raise ValueError("off_target needs guides — run crispr_design first.")
    seq = _require_sequence(step, ctx)
    references = [seq, *step.params.get("reference_sequences", [])]
    cas_family = step.params.get("cas_family", "cas9")
    max_mm = step.params.get("max_mismatches", 3)
    min_spec = step.params.get("min_specificity")  # optional pass/fail threshold

    annotated: list[dict] = []
    for g in guides:
        r = scan_off_targets(
            g["sequence"], references, cas_family=cas_family, max_mismatches=max_mm
        )
        g = dict(g)
        g["off_target_count"] = r.total
        g["off_target_summary"] = r.summary
        g["specificity_score"] = r.specificity_score
        annotated.append(g)

    passed = annotated
    if min_spec is not None:
        passed = [g for g in annotated if g["specificity_score"] >= min_spec]
    ctx["guides"] = passed
    return {"scanned": len(annotated), "passed": len(passed), "guides": passed}


def _h_primer_design(step: WorkflowStep, ctx: dict) -> dict:
    from ..routers.primers import PrimerDesignRequest, design_primers

    seq = _require_sequence(step, ctx)
    payload = dict(step.params)
    payload.pop("sequence", None)
    payload.setdefault("template", seq)
    pairs = design_primers(PrimerDesignRequest(**payload))
    ctx["primers"] = [p.model_dump() for p in pairs]
    return {"count": len(pairs), "primers": ctx["primers"]}


def _h_restriction_digest(step: WorkflowStep, ctx: dict) -> dict:
    from ..routers.primers import DigestRequest, digest

    seq = _require_sequence(step, ctx)
    enzymes = step.params.get("enzymes")
    if not enzymes:
        raise ValueError("restriction_digest needs params.enzymes (a list).")
    res = digest(DigestRequest(template=seq, enzymes=enzymes))
    ctx["digest"] = res.model_dump()
    return ctx["digest"]


def _h_export(step: WorkflowStep, ctx: dict) -> dict:
    source = step.params.get("source", "sequence")
    if source == "sequence":
        seq = _require_sequence(step, ctx)
        name = step.params.get("name", "sequence")
        content = f">{name}\n{seq}\n"
        filename, fmt = f"{name}.fasta", "fasta"
    elif source == "guides":
        guides = ctx.get("guides") or []
        header = "sequence\tpam\tposition\tstrand\tgc_content\ton_target_score\tspecificity\n"
        rows = "".join(
            f"{g.get('sequence', '')}\t{g.get('pam', '')}\t{g.get('position', '')}\t"
            f"{g.get('strand', '')}\t{g.get('gc_content', '')}\t"
            f"{g.get('on_target_score', '')}\t{g.get('specificity_score', '')}\n"
            for g in guides
        )
        content = header + rows
        filename, fmt = "guides.tsv", "tsv"
    elif source == "primers":
        primers = ctx.get("primers") or []
        header = "forward\treverse\tproduct_size\tforward_tm\treverse_tm\n"
        rows = "".join(
            f"{p.get('forward', {}).get('sequence', '')}\t"
            f"{p.get('reverse', {}).get('sequence', '')}\t{p.get('product_size', '')}\t"
            f"{p.get('forward', {}).get('tm', '')}\t{p.get('reverse', {}).get('tm', '')}\n"
            for p in primers
        )
        content = header + rows
        filename, fmt = "primers.tsv", "tsv"
    else:
        raise ValueError(f"export: unknown source '{source}' (sequence|guides|primers).")

    result = {"filename": filename, "format": fmt, "content": content}
    ctx["export"] = result
    return result


def _h_order(step: WorkflowStep, ctx: dict) -> dict:
    from ..formats import export_order
    from ..models.sequence import Sequence

    vendor = step.params.get("vendor")
    if not vendor:
        raise ValueError("order needs params.vendor (e.g. 'idt', 'twist', 'genewiz').")
    seq = _require_sequence(step, ctx)
    name = step.params.get("name", "order")
    start = step.params.get("start", 0)
    end = step.params.get("end")
    # Reuses the offline vendor file export — no order is submitted anywhere.
    data, filename, instructions = export_order(
        Sequence(id=name, name=name, seq=seq), vendor, start, end
    )
    result = {
        "vendor": vendor,
        "filename": filename,
        "file_bytes": len(data),
        "instructions": instructions,
        "submitted": False,
    }
    ctx["order"] = result
    return result


register_handler(StepType.CRISPR_DESIGN, _h_crispr_design)
register_handler(StepType.OFF_TARGET, _h_off_target)
register_handler(StepType.PRIMER_DESIGN, _h_primer_design)
register_handler(StepType.RESTRICTION_DIGEST, _h_restriction_digest)
register_handler(StepType.EXPORT, _h_export)
register_handler(StepType.ORDER, _h_order)


# ── Execution ────────────────────────────────────────────────────────────────

def run_step(step: WorkflowStep, ctx: dict) -> WorkflowStep:
    """Run a single step, updating its status/result/error in place."""
    step.status = StepStatus.RUNNING
    handler = HANDLERS.get(step.type)
    if handler is None:
        step.status = StepStatus.FAILED
        step.error = (
            f"Step type '{step.type.value}' is not supported by this engine yet "
            "(GUI/network steps are a follow-up)."
        )
        return step
    try:
        step.result = handler(step, ctx)
        step.status = StepStatus.COMPLETE
        step.error = None
    except HTTPException as e:  # operation-level validation error
        step.status = StepStatus.FAILED
        step.error = f"{e.status_code}: {e.detail}"
    except Exception as e:
        step.status = StepStatus.FAILED
        step.error = str(e)
    return step


def run_workflow(
    wf: Workflow, context: dict | None = None, *, stop_on_failure: bool = True
) -> dict:
    """Execute all steps sequentially, threading a shared context.

    Returns the final context. The workflow is mutated in place (each step's
    status/result/error). On failure, remaining steps are marked ``skipped``
    when ``stop_on_failure`` is set.
    """
    ctx: dict = dict(context or {})
    wf.reset()
    for i, step in enumerate(wf.steps):
        run_step(step, ctx)
        if step.status == StepStatus.FAILED and stop_on_failure:
            for later in wf.steps[i + 1:]:
                later.status = StepStatus.SKIPPED
            break
    return ctx
