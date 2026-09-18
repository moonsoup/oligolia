#!/usr/bin/env python3
"""
Phases 4+5: Compare & Locate
Diffs expected vs actual outputs, classifies failures, traces each one
back to a specific backend file and line number.

Input:  .claude/qa/runs/<latest>/actual_outputs.json
Output: .claude/qa/runs/<latest>/findings.json

Usage:
    python .claude/qa/analyst.py [run_id]
    python .claude/qa/analyst.py          # uses most recent run
"""

import json
import sys
import time
from pathlib import Path

BASE = Path(__file__).parent
RUNS_DIR = BASE / "runs"
BACKEND_DIR = Path(__file__).parent.parent.parent / "backend" / "routers"


def latest_run() -> Path:
    runs = sorted(RUNS_DIR.glob("*/actual_outputs.json"), key=lambda p: p.parent.name)
    if not runs:
        raise FileNotFoundError("No runs found. Run runner.py first.")
    return runs[-1]


# ── Route → source file map ───────────────────────────────────────────────────

def _discover_route_map() -> dict[str, str]:
    """Build the endpoint -> router-file map from the live app.

    It used to be a hand-maintained dict, and it had drifted: it named
    `/variants/blast`, which does not exist, and had no entries for `/cloning`
    or `/structure`. A map that is wrong in both directions makes the analyst
    point at the wrong file, or at no file at all (#69.3).

    Derived from `app.routes` so it cannot drift. The static table below is only
    a fallback for when the backend package cannot be imported.
    """
    try:
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
        from backend.main import app
    except Exception:
        return dict(_STATIC_ROUTE_MAP)

    routers_dir = Path(__file__).resolve().parents[2] / "backend" / "routers"
    discovered: dict[str, str] = {}
    for route in app.routes:
        path = getattr(route, "path", "")
        if path.count("/") < 2:
            continue
        prefix = path.lstrip("/").split("/")[0]
        module = f"{prefix}.py"
        # Only map paths that really belong to a router. FastAPI's own /docs and
        # /openapi.json, and any static mount, have no router file behind them.
        if not (routers_dir / module).is_file():
            continue
        discovered[path] = module
    return discovered or dict(_STATIC_ROUTE_MAP)


#: Fallback only. Kept so the analyst still works with no importable backend,
#: but never the source of truth — that is what let it go stale.
_STATIC_ROUTE_MAP = {
    "/sequences/":             "sequences.py",
    "/crispr/design":          "crispr.py",
    "/primers/design":         "primers.py",
    "/alignment/pairwise":     "alignment.py",
    "/alignment/multiple":     "alignment.py",
    "/variants/annotate":      "variants.py",
    "/analysis/composition":   "analysis.py",
    "/analysis/find_repeats":  "analysis.py",
    "/cloning/gibson":         "cloning.py",
    "/structure/predict":      "structure.py",
}

ROUTE_MAP = _discover_route_map()


def source_file_for(path: str) -> Path | None:
    for prefix, fname in ROUTE_MAP.items():
        if path.startswith(prefix):
            p = BACKEND_DIR / fname
            return p if p.exists() else None
    return None


def grep_source(path: Path, terms: list[str]) -> list[dict]:
    """Find lines in source that contain any search term."""
    if not path or not path.exists():
        return []
    hits = []
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines, 1):
        if any(t.lower() in line.lower() for t in terms):
            hits.append({"line": i, "content": line.rstrip()})
    return hits[:5]  # cap at 5 most relevant lines


# ── Validation helpers ────────────────────────────────────────────────────────

def check_result(actual: dict, expected: dict, test_id: str) -> list[dict]:
    """Return list of failures for a single test case."""
    failures = []
    body = actual.get("body")
    status = actual.get("status")
    error = actual.get("error")

    if error:
        failures.append({
            "type": "network_error",
            "detail": error,
            "severity": "error",
        })
        return failures

    # Status check
    exp_status = expected.get("status", 200)
    if isinstance(exp_status, list):
        if status not in exp_status:
            failures.append({
                "type": "wrong_status",
                "expected": exp_status,
                "actual": status,
                "severity": "critical" if status >= 500 else "major",
            })
    # A case may name several acceptable statuses. The MSA endpoint returns 200
    # with an aligner installed and an honest 503 without one (#58), and both are
    # correct -- reporting the 503 as a defect would make the pipeline argue
    # against a fix it should be confirming.
    elif status not in (expected.get("accept_status") or [exp_status]):
        failures.append({
            "type": "wrong_status",
            "expected": expected.get("accept_status") or exp_status,
            "actual": status,
            "severity": "critical" if status and status >= 500 else "major",
        })

    if body is None:
        return failures

    # Body-shape checks describe the SUCCESS schema. When a case accepts an
    # alternative status -- e.g. the MSA endpoint's honest 503 when no aligner is
    # installed (#58) -- the body is an error detail, so checking it for success
    # fields reports the correct behaviour as a missing field.
    if status != exp_status and status in (expected.get("accept_status") or []):
        return failures

    # Body shape checks
    if expected.get("is_list") and not isinstance(body, list):
        failures.append({
            "type": "wrong_type",
            "detail": f"Expected list, got {type(body).__name__}",
            "severity": "major",
        })

    if expected.get("body_contains") and isinstance(body, dict):
        for key in expected["body_contains"]:
            if key not in body:
                failures.append({
                    "type": "missing_key",
                    "key": key,
                    "severity": "major",
                })

    # Per-item checks on list responses
    if expected.get("each_item_has") and isinstance(body, list) and body:
        for field in expected["each_item_has"]:
            if field not in body[0]:
                failures.append({
                    "type": "missing_item_field",
                    "field": field,
                    "severity": "major",
                })

    # Numeric range checks
    if expected.get("length") and isinstance(body, dict):
        actual_len = body.get("length")
        if actual_len != expected["length"]:
            failures.append({
                "type": "wrong_value",
                "field": "length",
                "expected": expected["length"],
                "actual": actual_len,
                "severity": "minor",
            })

    if expected.get("gc_range") and isinstance(body, dict):
        gc = body.get("gc_content")
        if gc is not None:
            lo, hi = expected["gc_range"]
            if not (lo <= gc <= hi):
                failures.append({
                    "type": "value_out_of_range",
                    "field": "gc_content",
                    "value": gc,
                    "range": expected["gc_range"],
                    "severity": "major",
                })

    if expected.get("score_range") and isinstance(body, dict):
        guides = body.get("guides", [])
        for g in guides:
            score = g.get("score")
            if score is not None:
                lo, hi = expected["score_range"]
                if not (lo <= score <= hi):
                    failures.append({
                        "type": "value_out_of_range",
                        "field": "guide.score",
                        "value": score,
                        "range": expected["score_range"],
                        "severity": "major",
                    })
                    break

    if expected.get("guides_min") is not None and isinstance(body, dict):
        n = len(body.get("guides", []))
        if n < expected["guides_min"]:
            failures.append({
                "type": "too_few_results",
                "field": "guides",
                "expected_min": expected["guides_min"],
                "actual": n,
                "severity": "minor",
            })

    if expected.get("identity_range") and isinstance(body, dict):
        identity = body.get("identity")
        if identity is not None:
            lo, hi = expected["identity_range"]
            if not (lo <= identity <= hi):
                failures.append({
                    "type": "value_out_of_range",
                    "field": "identity",
                    "value": identity,
                    "range": expected["identity_range"],
                    "severity": "major",
                })

    # Performance
    latency = actual.get("latency_ms", 0)
    if latency > 5000:
        failures.append({
            "type": "performance",
            "latency_ms": latency,
            "threshold_ms": 5000,
            "severity": "warning",
        })

    return failures


SEVERITY_ORDER = {"critical": 0, "error": 1, "major": 2, "minor": 3, "warning": 4}


def locate_bug(path: str, failure: dict) -> dict:
    """Try to pinpoint the source line responsible for a failure."""
    src = source_file_for(path)
    if not src:
        return {"file": "unknown", "lines": []}

    search_terms = []
    ftype = failure.get("type", "")

    if ftype == "wrong_status":
        search_terms = ["raise", "HTTPException", "status_code"]
    elif ftype in ("missing_key", "missing_item_field"):
        key = failure.get("key") or failure.get("field", "")
        search_terms = [key, "response_model", "return"]
    elif ftype == "value_out_of_range":
        field = failure.get("field", "")
        search_terms = [field.split(".")[0], "score", "identity", "gc"]
    elif ftype == "too_few_results":
        search_terms = ["filter", "append", "guides", "result"]
    elif ftype == "wrong_type":
        search_terms = ["return", "list", "dict"]

    hits = grep_source(src, search_terms) if search_terms else []
    return {
        "file": str(src.relative_to(src.parent.parent.parent)),
        "lines": hits,
    }


def main():
    run_path = Path(sys.argv[1]) / "actual_outputs.json" if len(sys.argv) > 1 else latest_run()

    print("\n=== Phase 4+5: Compare & Locate ===")
    print(f"Analysing: {run_path}\n")

    with open(run_path) as f:
        run = json.load(f)

    findings = []

    for result in run["results"]:
        tc = result["test_case"]
        actual = result["actual"]
        expected = tc.get("expected", {})

        failures = check_result(actual, expected, tc["id"])
        if not failures:
            continue

        # Sort by severity
        failures.sort(key=lambda f: SEVERITY_ORDER.get(f.get("severity", "minor"), 99))
        worst = failures[0]["severity"]

        _, path = tc["endpoint"].split(" ", 1)
        location = locate_bug(path, failures[0])

        finding = {
            "test_id": tc["id"],
            "endpoint": tc["endpoint"],
            "tags": tc.get("tags", []),
            "known_bug": tc.get("expected", {}).get("known_bug"),
            "severity": worst,
            "failures": failures,
            "source_location": location,
            "actual_status": actual.get("status"),
            "actual_body_preview": str(actual.get("body", ""))[:300],
            "latency_ms": actual.get("latency_ms"),
        }
        findings.append(finding)

        icon = {"critical": "💥", "error": "✗", "major": "✗", "minor": "!", "warning": "⚠"}.get(worst, "?")
        print(f"  {icon} [{worst:8}] {tc['id']}")
        for f in failures[:2]:
            print(f"             {f['type']}: {f.get('detail') or f.get('expected','') or f.get('field','')}")
        if location["lines"]:
            print(f"             → {location['file']}:{location['lines'][0]['line']}")

    out = run_path.parent / "findings.json"
    with open(out, "w") as f:
        json.dump({
            "run_id": run["run_id"],
            "analysed": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "total_failures": len(findings),
            "by_severity": {
                sev: len([x for x in findings if x["severity"] == sev])
                for sev in ["critical", "error", "major", "minor", "warning"]
            },
            "findings": sorted(findings, key=lambda x: SEVERITY_ORDER.get(x["severity"], 99)),
        }, f, indent=2)

    print(f"\n{'='*55}")
    print(f"  {len(findings)} findings → {out}")
    print("  Severity: " + "  ".join(
        f"{sev}={len([x for x in findings if x['severity']==sev])}"
        for sev in ["critical", "major", "minor", "warning"]
    ))

    # Update agent_comms
    comms_path = BASE.parent / "agent_comms.json"
    if comms_path.exists():
        with open(comms_path) as f:
            obj = json.load(f)
        obj["workflow_state"]["current_phase"] = 6
        for phase in ["compare", "locate"]:
            if phase not in obj["workflow_state"]["phases_completed"]:
                obj["workflow_state"]["phases_completed"].append(phase)
        for msg in obj["messages"]:
            if msg.get("id") == "msg_2":
                msg["status"] = "accepted"
        obj["messages"].append({
            "id": "msg_3",
            "from": "analyst",
            "to": "reporter",
            "type": "handoff",
            "subject": f"{len(findings)} findings ready to report",
            "status": "pending",
            "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "body": (
                f"Analysis complete. {len(findings)} findings at "
                f".claude/qa/runs/{run['run_id']}/findings.json. "
                f"Run reporter.py to file GitHub issues for unfiled findings."
            ),
        })
        with open(comms_path, "w") as f:
            json.dump(obj, f, indent=2)
        print("✓ agent_comms.json updated")


if __name__ == "__main__":
    main()
