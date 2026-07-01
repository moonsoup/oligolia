"""Tests for the multi-step workflow engine (backend/workflow/)."""

from ..workflow import (
    StepStatus,
    StepType,
    Workflow,
    WorkflowStep,
    dumps,
    load_ogo,
    loads,
    run_step,
    run_workflow,
    save_ogo,
)


def _wf(*steps: WorkflowStep) -> Workflow:
    return Workflow(name="test", steps=list(steps))


def test_crispr_then_offtarget_pipeline(tp53_exon7: str) -> None:
    """Sequence feeds CRISPR design, whose guides feed the off-target step."""
    wf = _wf(
        WorkflowStep(id="s1", type=StepType.CRISPR_DESIGN,
                     params={"sequence": tp53_exon7, "max_guides": 5}),
        WorkflowStep(id="s2", type=StepType.OFF_TARGET, params={"min_specificity": 0}),
    )
    ctx = run_workflow(wf)
    assert wf.step("s1").status == StepStatus.COMPLETE
    assert wf.step("s2").status == StepStatus.COMPLETE
    assert wf.is_complete
    # Off-target step annotated every guide with a specificity score.
    guides = ctx["guides"]
    assert guides and all("specificity_score" in g for g in guides)


def test_context_threads_sequence_forward(tp53_exon7: str) -> None:
    """A sequence produced/held in context reaches a later primer step."""
    wf = _wf(
        WorkflowStep(id="s1", type=StepType.CRISPR_DESIGN,
                     params={"sequence": tp53_exon7, "max_guides": 3}),
        WorkflowStep(id="s2", type=StepType.PRIMER_DESIGN, params={"product_min": 60}),
    )
    ctx = run_workflow(wf)
    assert wf.step("s2").status == StepStatus.COMPLETE
    assert "primers" in ctx


def test_restriction_digest_step(tp53_exon7: str) -> None:
    wf = _wf(
        WorkflowStep(id="s1", type=StepType.RESTRICTION_DIGEST,
                     params={"sequence": tp53_exon7, "enzymes": ["EcoRI", "BamHI"]}),
    )
    run_workflow(wf)
    step = wf.step("s1")
    assert step.status == StepStatus.COMPLETE
    assert step.result["template_length"] == len(tp53_exon7)


def test_export_guides_step(tp53_exon7: str) -> None:
    wf = _wf(
        WorkflowStep(id="s1", type=StepType.CRISPR_DESIGN,
                     params={"sequence": tp53_exon7, "max_guides": 4}),
        WorkflowStep(id="s2", type=StepType.EXPORT, params={"source": "guides"}),
    )
    ctx = run_workflow(wf)
    exp = ctx["export"]
    assert exp["filename"] == "guides.tsv"
    assert exp["content"].startswith("sequence\tpam\t")
    assert len(exp["content"].splitlines()) >= 2  # header + at least one guide


def test_order_step_exports_no_submission(tp53_exon7: str) -> None:
    """The order step produces a vendor file and never submits an order."""
    wf = _wf(
        WorkflowStep(id="s1", type=StepType.ORDER,
                     params={"sequence": tp53_exon7[:120], "vendor": "twist", "name": "myorder"}),
    )
    run_workflow(wf)
    step = wf.step("s1")
    assert step.status == StepStatus.COMPLETE
    assert step.result["submitted"] is False
    assert step.result["file_bytes"] > 0


def test_failure_stops_and_skips_remaining(tp53_exon7: str) -> None:
    """A failing step halts the run and marks later steps skipped."""
    wf = _wf(
        WorkflowStep(id="s1", type=StepType.RESTRICTION_DIGEST,
                     params={"sequence": tp53_exon7, "enzymes": ["NotAnEnzyme"]}),
        WorkflowStep(id="s2", type=StepType.EXPORT, params={"source": "sequence"}),
    )
    run_workflow(wf)
    assert wf.step("s1").status == StepStatus.FAILED
    assert wf.step("s1").error
    assert wf.step("s2").status == StepStatus.SKIPPED


def test_missing_sequence_reports_clear_error() -> None:
    wf = _wf(WorkflowStep(id="s1", type=StepType.PRIMER_DESIGN, params={}))
    run_workflow(wf)
    step = wf.step("s1")
    assert step.status == StepStatus.FAILED
    assert "No sequence available" in step.error


def test_unsupported_step_type_is_non_fatal(tp53_exon7: str) -> None:
    """db_search isn't wired in this slice; it fails clearly rather than crashing."""
    wf = _wf(WorkflowStep(id="s1", type=StepType.DB_SEARCH, params={"query": "TP53"}))
    run_workflow(wf)
    step = wf.step("s1")
    assert step.status == StepStatus.FAILED
    assert "not supported" in step.error


def test_ogo_roundtrip(tmp_path, tp53_exon7: str) -> None:
    """A run workflow serializes to .ogo and reloads with results intact."""
    wf = _wf(
        WorkflowStep(id="s1", type=StepType.CRISPR_DESIGN,
                     params={"sequence": tp53_exon7, "max_guides": 3}),
    )
    run_workflow(wf)
    path = tmp_path / "design.ogo"
    save_ogo(wf, path)
    reloaded = load_ogo(path)
    assert reloaded.name == wf.name
    assert reloaded.step("s1").status == StepStatus.COMPLETE
    assert reloaded.step("s1").result["total_candidates"] == wf.step("s1").result["total_candidates"]
    # String round-trip too.
    assert loads(dumps(wf)).step("s1").status == StepStatus.COMPLETE


def test_reset_clears_run_state(tp53_exon7: str) -> None:
    wf = _wf(
        WorkflowStep(id="s1", type=StepType.CRISPR_DESIGN,
                     params={"sequence": tp53_exon7, "max_guides": 2}),
    )
    run_workflow(wf)
    assert wf.step("s1").status == StepStatus.COMPLETE
    wf.reset()
    assert wf.step("s1").status == StepStatus.PENDING
    assert wf.step("s1").result is None


def test_run_single_step_directly(tp53_exon7: str) -> None:
    step = WorkflowStep(id="s1", type=StepType.CRISPR_DESIGN,
                        params={"sequence": tp53_exon7, "max_guides": 2})
    ctx: dict = {}
    run_step(step, ctx)
    assert step.status == StepStatus.COMPLETE
    assert ctx["guides"]
