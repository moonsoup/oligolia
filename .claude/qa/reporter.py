#!/usr/bin/env python3
"""
Phase 6: Report
Reads findings.json, deduplicates against open GitHub issues,
files new issues for unfiled bugs via gh CLI.

Usage:
    python .claude/qa/reporter.py [run_id]
    python .claude/qa/reporter.py --dry-run
"""

import hashlib
import os
import tempfile
import json
import sys
import time
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import _state  # noqa: E402

BASE = Path(__file__).parent
RUNS_DIR = BASE / "runs"
REPO = "moonsoup/oligolia"

DRY_RUN = "--dry-run" in sys.argv


def latest_findings() -> Path:
    runs = sorted(RUNS_DIR.glob("*/findings.json"), key=lambda p: p.parent.name)
    if not runs:
        raise FileNotFoundError("No findings. Run analyst.py first.")
    return runs[-1]


def open_issues() -> list[dict]:
    """Fetch open GitHub issues as list of {number, title, body}."""
    result = subprocess.run(
        ["gh", "issue", "list", "--repo", REPO, "--state", "open",
         "--limit", "100", "--json", "number,title,body"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"  WARN: Could not fetch open issues: {result.stderr}")
        return []
    return json.loads(result.stdout)


def fingerprint(finding: dict) -> str:
    """A stable identity for a finding: endpoint + defect class + assertion.

    The old rule was "two or more shared lowercase words with an open issue
    title". With issue titles full of shared vocabulary -- crispr, guide,
    sequence, primer, variant, alignment -- a genuinely new finding on a
    different endpoint collided with an unrelated issue and was silently
    dropped. A discovery pipeline that discards findings without saying so has
    the same shape as the fallback it exists to catch (#74).
    """
    parts = [
        finding.get("endpoint", ""),
        finding.get("category", finding.get("severity", "")),
        finding.get("assertion", finding.get("title", "")),
    ]
    digest = hashlib.sha256("|".join(str(p) for p in parts).encode()).hexdigest()[:12]
    return f"qa-fingerprint:{digest}"


def is_duplicate(finding: dict, open_issues: list[dict]) -> bool:
    """Has this exact finding already been filed?

    Matched on the fingerprint recorded in the issue body, so it is an identity
    check rather than a guess from overlapping words.
    """
    fp = fingerprint(finding)
    for issue in open_issues:
        if fp in (issue.get("body") or ""):
            return True
    return False


def format_issue(finding: dict) -> tuple[str, str]:
    """Format a finding as a GitHub issue title + body."""
    endpoint = finding["endpoint"]
    failures = finding["failures"]
    loc = finding["source_location"]
    sev = finding["severity"]
    tags = finding.get("tags", [])

    primary = failures[0]
    ftype = primary.get("type", "unknown")

    # Build title
    if ftype == "wrong_status":
        title = f"[{sev.upper()}] {endpoint} returns {finding['actual_status']} instead of {primary['expected']}"
    elif ftype == "missing_key":
        title = f"[{sev.upper()}] {endpoint} response missing field `{primary.get('key')}`"
    elif ftype == "value_out_of_range":
        title = f"[{sev.upper()}] {endpoint} `{primary.get('field')}` value out of expected range"
    elif ftype == "performance":
        title = f"[PERF] {endpoint} response time {finding['latency_ms']}ms exceeds 5s threshold"
    elif ftype == "network_error":
        title = f"[{sev.upper()}] {endpoint} causes unhandled exception"
    else:
        title = f"[{sev.upper()}] {endpoint} — {ftype.replace('_', ' ')}"

    # Build body
    lines = [
        f"## Auto-detected by QA runner (run `{finding.get('test_id','')}` in `.claude/qa/`)\n",
        f"**Endpoint:** `{endpoint}`  ",
        f"**Severity:** {sev}  ",
        f"**Tags:** {', '.join(tags)}\n",
        "## Failures\n",
    ]
    for f in failures:
        detail = " | ".join(f"{k}: {v}" for k, v in f.items() if k != "type" and v is not None)
        lines.append(f"- **{f['type']}**: {detail}")

    lines += [
        "\n## Observed Response\n",
        f"```\nHTTP {finding['actual_status']}  ({finding['latency_ms']}ms)\n",
        finding.get("actual_body_preview", "(no body)"),
        "```\n",
    ]

    if loc.get("file") and loc.get("lines"):
        lines += ["\n## Likely Source\n", f"`{loc['file']}`\n"]
        for hit in loc["lines"][:3]:
            lines.append(f"- Line {hit['line']}: `{hit['content'][:100]}`")

    if finding.get("known_bug"):
        lines += [f"\n> **Note:** This matches known bug `{finding['known_bug']}`"]

    # The dedup identity, recorded in the body so a later run matches on it
    # rather than guessing from overlapping title words (#74).
    lines += ["", f"<!-- {fingerprint(finding)} -->"]
    return title, "\n".join(lines)


#: The project's disciplined filing interface. It refuses an empty body, takes
#: the body only from a file rather than an inline string, and refuses to close
#: an issue without evidence or one assigned to someone else. Calling
#: `gh issue create` directly bypassed all of that (#74).
GH_ISSUE = Path.home() / "Software" / "projectMan" / "scripts" / "gh_issue.py"


def file_issue(title: str, body: str, labels: list[str]) -> str | None:
    if DRY_RUN:
        print(f"\n  [DRY RUN] Would file: {title}")
        return None

    if not GH_ISSUE.is_file():
        print(f"  WARN: {GH_ISSUE} not found; refusing to file directly with gh (#74)")
        return None

    # gh_issue.py takes the body from a file, never an inline argument.
    with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False) as fh:
        fh.write(body)
        body_path = fh.name

    try:
        result = subprocess.run(
            [sys.executable, str(GH_ISSUE), "file",
             "--repo", REPO, "--title", title, "--body-file", body_path],
            capture_output=True, text=True,
        )
    finally:
        os.unlink(body_path)

    # 0 = filed, 2 = refused by its own rules, 3 = gh failed.
    if result.returncode == 0:
        for line in result.stdout.splitlines():
            if line.startswith("FILED:"):
                return line.split(":", 1)[1].strip()
        return result.stdout.strip()
    if result.returncode == 2:
        print(f"  REFUSED by gh_issue.py: {result.stdout.strip() or result.stderr.strip()}")
    else:
        print(f"  WARN: filing failed ({result.returncode}): {result.stderr.strip()}")
    return None


def main():
    findings_path = latest_findings()
    print("\n=== Phase 6: Report ===")
    print(f"Findings: {findings_path}\n")

    with open(findings_path) as f:
        data = json.load(f)

    findings = data["findings"]
    existing = open_issues()
    print(f"  {len(existing)} open issues fetched for dedup check\n")

    filed = []
    skipped = []

    for finding in findings:
        sev = finding["severity"]

        # Skip warnings and minor from auto-filing (keep signal high)
        if sev in ("warning", "minor"):
            skipped.append((finding["test_id"], "low severity"))
            continue

        # Skip known bugs that are already filed
        if finding.get("known_bug"):
            skipped.append((finding["test_id"], "known bug already tracked"))
            continue

        if is_duplicate(finding, existing):
            skipped.append((finding["test_id"], "likely duplicate"))
            continue

        title, body = format_issue(finding)
        labels = ["bug"] if sev in ("critical", "error", "major") else ["enhancement"]

        url = file_issue(title, body, labels)
        if url:
            filed.append((finding["test_id"], url))
            print(f"  ✓ Filed: {url}")
            time.sleep(1)  # GH rate limit
        else:
            if DRY_RUN:
                filed.append((finding["test_id"], "[dry-run]"))

    print(f"\n{'='*55}")
    print(f"  Filed: {len(filed)}  Skipped: {len(skipped)}")
    for tid, reason in skipped:
        print(f"    skip  {tid}: {reason}")

    # Update agent_comms
    obj = _state.load()
    # The last phase. `None` meant "finished" but read as "unknown" in the
    # state file; len(PHASES) is unambiguous.
    obj["workflow_state"]["current_phase"] = len(_state.PHASES)
    obj["workflow_state"]["last_run"] = data["run_id"]
    for phase in ["report"]:
        if phase not in obj["workflow_state"]["phases_completed"]:
            obj["workflow_state"]["phases_completed"].append(phase)
    for msg in obj["messages"]:
        if msg.get("id") == "msg_3":
            msg["status"] = "accepted"
    obj["messages"].append({
        "id": "msg_4",
        "from": "reporter",
        "to": "coordinator",
        "type": "status_update",
        "subject": f"Run {data['run_id']} complete — {len(filed)} issues filed",
        "status": "resolved",
        "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "body": (
            f"All 6 phases complete. Filed {len(filed)} new issues. "
            f"Skipped {len(skipped)} (duplicates, low severity, known bugs). "
            f"Full findings: .claude/qa/runs/{data['run_id']}/findings.json"
        ),
    })
    _state.save(obj)
    print("✓ pipeline_state.json updated — all phases complete")


if __name__ == "__main__":
    main()
