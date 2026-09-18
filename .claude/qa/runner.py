#!/usr/bin/env python3
"""
Phase 3: Execution
Drives every test case in the corpus against the running Oligolia backend.
Records actual HTTP status, response body, latency, and stderr.

Requires: backend running on localhost:8765 (override with OLIGOLIA_PORT)
Output:   .claude/qa/runs/<timestamp>/actual_outputs.json

Usage:
    # Terminal 1: start backend
    cd /path/to/oligolia && python run_backend.py

    # Terminal 2: run tests
    python .claude/qa/runner.py
"""

import json
import os
import time
import urllib.request
import urllib.error
import urllib.parse
from pathlib import Path
from datetime import datetime

BASE = Path(__file__).parent
CORPUS_PATH = BASE / "corpus" / "corpus.json"
RUNS_DIR = BASE / "runs"
# run_backend.py serves 8765. This said 8000, so the runner could not connect
# even when the pipeline was present (#69.2). backend/tests/test_qa_pipeline.py
# asserts the two agree, so the next change to either is caught.
BACKEND_PORT = int(os.environ.get("OLIGOLIA_PORT", "8765"))
BACKEND_URL = f"http://localhost:{BACKEND_PORT}"


def wait_for_backend(timeout=30):
    print(f"Waiting for backend at {BACKEND_URL}...")
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            urllib.request.urlopen(f"{BACKEND_URL}/health", timeout=2)
            print("  ✓ Backend is up\n")
            return True
        except Exception:
            time.sleep(1)
    # Was a plain string, so it printed the literal "{timeout}s".
    print(f"  ✗ Backend not reachable after {timeout}s — is run_backend.py running?")
    return False


def call_endpoint(test_case: dict) -> dict:
    endpoint = test_case["endpoint"]
    method, path = endpoint.split(" ", 1)
    url = BACKEND_URL + path

    # Several endpoints take their input as QUERY parameters rather than a JSON
    # body -- /analysis/composition, find_orfs and protein_properties all declare
    # `sequence: str` as a function argument, which FastAPI reads from the query
    # string. The corpus posted them as a body and got 422 "Field required"
    # against ('query', 'sequence') every time (#69).
    if test_case.get("query"):
        url += "?" + urllib.parse.urlencode(test_case["query"])

    start = time.time()
    result = {
        "id": test_case["id"],
        "endpoint": endpoint,
        "status": None,
        "body": None,
        "error": None,
        "latency_ms": None,
        "ts": datetime.utcnow().isoformat() + "Z",
    }

    try:
        # VCF uploads use multipart — simplified to raw text upload for now
        if test_case.get("payload_type") == "vcf_upload":
            boundary = "----OligoliaBoundary"
            body = (
                f"--{boundary}\r\n"
                f'Content-Disposition: form-data; name="file"; filename="test.vcf"\r\n'
                f"Content-Type: text/plain\r\n\r\n"
                + test_case["payload_raw"]
                + f"\r\n--{boundary}--\r\n"
            ).encode("utf-8")
            req = urllib.request.Request(
                url, data=body, method="POST",
                headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
            )
        elif method == "POST":
            payload = json.dumps(test_case.get("payload", {})).encode()
            req = urllib.request.Request(
                url, data=payload, method="POST",
                headers={"Content-Type": "application/json"},
            )
        else:
            req = urllib.request.Request(url, method="GET")

        with urllib.request.urlopen(req, timeout=30) as r:
            result["status"] = r.status
            result["body"] = json.loads(r.read())

    except urllib.error.HTTPError as e:
        result["status"] = e.code
        try:
            result["body"] = json.loads(e.read())
        except Exception:
            result["body"] = None
    except Exception as e:
        result["error"] = str(e)

    result["latency_ms"] = round((time.time() - start) * 1000)
    return result


def main():
    # `return` here exited 0, so run_all.sh's `set -e` could not stop the
    # pipeline: the analyst then ran with no run directory and died with a
    # traceback instead of a message (#69).
    if not CORPUS_PATH.exists():
        print(f"✗ Corpus not found at {CORPUS_PATH}. Run scout.py first.")
        raise SystemExit(2)

    with open(CORPUS_PATH) as f:
        corpus = json.load(f)

    if not wait_for_backend():
        raise SystemExit(3)

    run_id = datetime.utcnow().strftime("%Y%m%dT%H%M%S")
    run_dir = RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    test_cases = corpus["test_cases"]
    print(f"Running {len(test_cases)} test cases...\n")

    results = []
    passed = failed = errored = slow = 0

    for tc in test_cases:
        actual = call_endpoint(tc)
        expected = tc.get("expected", {})

        # Quick pass/fail for reporting
        exp_status = expected.get("status", 200)
        if isinstance(exp_status, list):
            ok = actual["status"] in exp_status
        else:
            ok = actual["status"] == exp_status

        if actual["error"]:
            status_str = "ERROR"
            errored += 1
        elif ok:
            status_str = "PASS"
            passed += 1
        else:
            status_str = "FAIL"
            failed += 1

        if actual["latency_ms"] and actual["latency_ms"] > 5000:
            slow += 1
            status_str += " (SLOW)"

        print(f"  [{status_str:12}] {tc['id']}  ({actual['latency_ms']}ms)  HTTP {actual['status']}")

        results.append({
            "test_case": tc,
            "actual": actual,
            "pass": ok and not actual["error"],
        })

        time.sleep(0.05)  # gentle pacing

    out = run_dir / "actual_outputs.json"
    with open(out, "w") as f:
        json.dump({
            "run_id": run_id,
            "corpus_version": corpus.get("version"),
            "backend_url": BACKEND_URL,
            "summary": {
                "total": len(test_cases),
                "passed": passed,
                "failed": failed,
                "errored": errored,
                "slow": slow,
            },
            "results": results,
        }, f, indent=2)

    print(f"\n{'='*55}")
    print(f"  Total: {len(test_cases)}  Pass: {passed}  Fail: {failed}  Error: {errored}  Slow: {slow}")
    print(f"  Results → {out}")

    # Update agent_comms
    comms_path = BASE.parent / "agent_comms.json"
    if comms_path.exists():
        with open(comms_path) as f:
            obj = json.load(f)
        obj["workflow_state"]["current_phase"] = 4
        if "execution" not in obj["workflow_state"]["phases_completed"]:
            obj["workflow_state"]["phases_completed"].append("execution")
        # Mark scout handoff accepted
        for msg in obj["messages"]:
            if msg.get("id") == "msg_1":
                msg["status"] = "accepted"
        obj["messages"].append({
            "id": "msg_2",
            "from": "runner",
            "to": "analyst",
            "type": "handoff",
            "subject": f"Execution complete — {failed} failures, {errored} errors",
            "status": "pending",
            "ts": datetime.utcnow().isoformat() + "Z",
            "body": (
                f"Run {run_id} complete. {passed} passed, {failed} failed, {errored} errors, {slow} slow. "
                f"Results at .claude/qa/runs/{run_id}/actual_outputs.json. "
                f"Run analyst.py to diff expected vs actual and locate bugs."
            ),
        })
        with open(comms_path, "w") as f:
            json.dump(obj, f, indent=2)
        print("✓ agent_comms.json updated")

    return run_id


if __name__ == "__main__":
    main()
