#!/usr/bin/env bash
# Full QA pipeline: scout → runner → analyst → reporter
# Run from the oligolia project root.
#
# Usage:
#   bash .claude/qa/run_all.sh           # full run, files issues
#   bash .claude/qa/run_all.sh --dry-run # full run, no issue filing
#
# -e stops on a failing phase, but a phase has to actually fail: runner.py used
# to `return` (exit 0) when the backend was unreachable, so the analyst ran with
# no input and died with a traceback (#69). The exit codes are real now.
set -euo pipefail
QA=".claude/qa"
DRY=${1:-""}

fail() { echo ""; echo "✗ $1"; exit "$2"; }

echo "=== Oligolia QA Pipeline ==="
echo ""

echo "▶ Phase 1+2: Scout (discovery + corpus generation)"
python "$QA/scout.py"

echo ""
echo "▶ Phase 3: Runner (execution against live backend)"
echo "   Make sure backend is running: python run_backend.py"
python "$QA/runner.py" || fail "runner failed (exit $?) — is the backend up on 8765?" 3

echo ""
echo "▶ Phase 4+5: Analyst (compare + locate)"
python "$QA/analyst.py"

echo ""
echo "▶ Phase 6: Reporter (file GitHub issues)"
python "$QA/reporter.py" $DRY

echo ""
echo "=== Done. Check .claude/qa/runs/ for full output. ==="
