"""#69: the documented QA pipeline could not run. Pinned so it cannot drift again.

Three defects, all confirmed:

  1. `.gitignore` ignored `.claude/qa/`, so the remote dev loop cloned main and
     got a repo with no pipeline — while CLAUDE.md told it to run
     `bash .claude/qa/run_all.sh` as its DISCOVER step;
  2. the runner polled port 8000; `run_backend.py` serves 8765;
  3. the endpoint map was hand-maintained and had drifted — it mapped
     `/variants/blast`, which does not exist, and had no `/cloning` or
     `/structure` entries.

`--dry-run` worked, which is why none of it was noticed. These tests are the
mechanical check that replaces noticing.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

QA = Path(__file__).resolve().parents[2] / ".claude" / "qa"

pytestmark = pytest.mark.skipif(
    not QA.is_dir(), reason=".claude/qa is not present in this checkout"
)


def _module_constant(path: Path, name: str):
    """Read a module-level literal without importing."""
    tree = ast.parse(path.read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == name:
                    return ast.literal_eval(node.value)
    raise AssertionError(f"{name} not found in {path.name}")


def _load(name: str):
    """Import a QA script by path. All four guard their entry point with
    `if __name__ == "__main__"`, so importing them is side-effect free."""
    import importlib.util

    import sys
    if str(QA) not in sys.path:
        sys.path.insert(0, str(QA))
    spec = importlib.util.spec_from_file_location(f"qa_{name}", QA / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_pipeline_is_tracked_in_git() -> None:
    """#69.1: the remote loop is told to run this, so it has to be in the repo."""
    gitignore = (QA.parents[1] / ".gitignore").read_text()
    offending = [
        line for line in gitignore.splitlines()
        if line.strip() and not line.strip().startswith("#")
        and line.strip().rstrip("/") in (".claude/qa", ".claude/qa/*")
    ]
    assert not offending, f".gitignore still excludes the QA pipeline: {offending}"


def test_the_runner_targets_the_port_the_backend_serves() -> None:
    """#69.2: the runner polled 8000 and the backend serves 8765."""
    backend_py = (QA.parents[1] / "run_backend.py").read_text()
    served = re.search(r"port\s*=\s*(\d+)", backend_py)
    assert served, "could not find the port in run_backend.py"

    url = _load("runner").BACKEND_URL
    assert served.group(1) in url, (
        f"runner targets {url} but run_backend.py serves port {served.group(1)}"
    )


def test_every_mapped_route_exists_in_the_app() -> None:
    """#69.3: the map listed /variants/blast, which does not exist."""
    from backend.main import app

    real = {getattr(r, "path", None) for r in app.routes}
    mapped = _load("analyst").ROUTE_MAP

    missing = sorted(p for p in mapped if p not in real)
    assert not missing, f"the endpoint map names routes the app does not serve: {missing}"


def test_the_router_modules_named_in_the_map_exist() -> None:
    mapped = _load("analyst").ROUTE_MAP
    routers = QA.parents[1] / "backend" / "routers"
    missing = sorted({m for m in mapped.values() if not (routers / m).is_file()})
    assert not missing, f"the map names router files that do not exist: {missing}"


def test_the_map_covers_the_routers_that_have_endpoints() -> None:
    """Drift the other way: a whole router missing from the map is invisible."""
    from backend.main import app

    served_prefixes = {
        "/" + getattr(r, "path", "/").lstrip("/").split("/")[0]
        for r in app.routes
        if getattr(r, "path", "/").count("/") >= 2
    }
    mapped_prefixes = {
        "/" + p.lstrip("/").split("/")[0]
        for p in _load("analyst").ROUTE_MAP
    }
    uncovered = sorted(served_prefixes - mapped_prefixes - {"/openapi.json", "/docs", "/redoc"})
    assert not uncovered, (
        f"these router prefixes have endpoints but no entry in the map: {uncovered}"
    )


def test_the_scout_does_not_substitute_a_stub_on_fetch_failure() -> None:
    """#73: a failed fetch became 92 nt of HBB under the requested gene's name."""
    source = (QA / "scout.py").read_text()
    assert "FALLBACK" not in source.upper() or "raise" in source, (
        "scout.py still has a fallback path; a failed fetch must skip or fail, "
        "not substitute another gene's sequence under the requested name (#73)"
    )


def test_the_reporter_files_through_the_disciplined_interface() -> None:
    """#74: it called `gh issue create` directly, bypassing the refusals."""
    source = (QA / "reporter.py").read_text()
    assert "gh_issue.py" in source, (
        "reporter.py must file through projectMan/scripts/gh_issue.py, which "
        "refuses an empty body and a close without evidence (#74)"
    )


def test_the_discovered_map_is_not_empty_and_beats_the_fallback() -> None:
    """The map must come from the app, not from the static fallback.

    If discovery silently failed, the fallback would keep the analyst working
    while going stale again — which is how this got here.
    """
    analyst = _load("analyst")
    static = analyst._STATIC_ROUTE_MAP
    assert analyst.ROUTE_MAP, "the route map is empty"
    assert len(analyst.ROUTE_MAP) > len(static), (
        f"discovery produced {len(analyst.ROUTE_MAP)} routes, no more than the "
        f"{len(static)}-entry fallback — discovery is probably not running"
    )


# --- #74: the dedup identity must be an identity, not a word overlap ---

def _finding(**over) -> dict:
    base = {
        "endpoint": "POST /crispr/design", "severity": "error",
        "category": "wrong_value", "assertion": "guide excludes PAM",
        "tags": ["crispr"],
        "failures": [{"type": "wrong_value", "detail": "guide included a PAM base"}],
        "source_location": {}, "actual_status": 200, "latency_ms": 5,
    }
    base.update(over)
    return base


def test_the_fingerprint_is_stamped_into_the_issue_body() -> None:
    reporter = _load("reporter")
    _title, body = reporter.format_issue(_finding())
    assert reporter.fingerprint(_finding()) in body


def test_a_finding_matches_its_own_previous_issue() -> None:
    reporter = _load("reporter")
    _title, body = reporter.format_issue(_finding())
    assert reporter.is_duplicate(_finding(), [{"number": 1, "title": "x", "body": body}])


def test_a_word_overlapping_but_unrelated_issue_is_not_a_duplicate() -> None:
    """#74: "two or more shared lowercase title words" suppressed real findings.

    With issue titles full of shared vocabulary — crispr, guide, sequence,
    primer, variant, alignment — a genuinely new finding collided with an
    unrelated issue and was dropped without a word.
    """
    reporter = _load("reporter")
    unrelated = {"number": 2, "title": "crispr guide sequence design is wrong", "body": "no marker"}
    assert not reporter.is_duplicate(_finding(), [unrelated])


def test_a_different_endpoint_is_a_different_finding() -> None:
    reporter = _load("reporter")
    _title, body = reporter.format_issue(_finding())
    other = _finding(endpoint="POST /primers/design")
    assert not reporter.is_duplicate(other, [{"number": 1, "title": "x", "body": body}])


def test_the_fingerprint_is_stable_across_runs() -> None:
    reporter = _load("reporter")
    assert reporter.fingerprint(_finding()) == reporter.fingerprint(_finding())


# --- #80: the pipeline's own state, defined once and actually written ---

def test_the_state_path_is_defined_in_one_place() -> None:
    """Ziggurat flagged agent_comms.json as one path with four namers."""
    for name in ("scout", "runner", "analyst", "reporter"):
        source = (QA / f"{name}.py").read_text()
        assert 'BASE.parent / "agent_comms.json"' not in source, (
            f"{name}.py still builds its own path to a file it does not own (#80)"
        )
        assert "_state" in source, f"{name}.py should use the shared state module"


def test_the_state_file_is_not_the_codex_protocol_file() -> None:
    """Two unrelated protocols must not share a filename.

    The repo root's agent_comms.json is the Claude<->Codex review protocol
    (`_protocol` + typed messages). The pipeline writes `workflow_state`. Pointing
    one at the other crashes on a missing key (#80).
    """
    state = _load("_state")
    assert state.STATE_PATH.name != "agent_comms.json"
    assert state.STATE_PATH.parent == QA, state.STATE_PATH


def test_the_state_round_trips(tmp_path, monkeypatch) -> None:
    state = _load("_state")
    monkeypatch.setattr(state, "STATE_PATH", tmp_path / "pipeline_state.json")

    fresh = state.load()
    assert fresh["workflow_state"]["current_phase"] == 0
    assert fresh["workflow_state"]["phases_completed"] == []

    state.record(3, ["discover", "corpus"], {"id": "msg_1", "from": "scout"})
    again = state.load()
    assert again["workflow_state"]["current_phase"] == 3
    assert again["workflow_state"]["phases_completed"] == ["discover", "corpus"]
    assert len(again["messages"]) == 1


def test_recording_a_phase_twice_does_not_duplicate_it(tmp_path, monkeypatch) -> None:
    state = _load("_state")
    monkeypatch.setattr(state, "STATE_PATH", tmp_path / "pipeline_state.json")
    state.record(1, ["discover"])
    state.record(2, ["discover", "corpus"])
    assert state.load()["workflow_state"]["phases_completed"] == ["discover", "corpus"]


def test_a_corrupt_state_file_is_replaced_not_fatal(tmp_path, monkeypatch) -> None:
    state = _load("_state")
    path = tmp_path / "pipeline_state.json"
    path.write_text("{not json")
    monkeypatch.setattr(state, "STATE_PATH", path)
    assert state.load()["workflow_state"]["current_phase"] == 0
