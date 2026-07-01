---
name: advance-queue
description: Run one full iteration of the Oligolia dev loop — orient on open issues, verify recently pushed fixes, discover new issues via the QA pipeline, health-check the remote.
user-invocable: true
---

# Oligolia — Dev Loop (one iteration)

Adapted from atp-siege's `/advance-queue` (see `~/Software/bio_battle/.claude/skills/advance-queue/SKILL.md`).
The remote Docker instance fixes GitHub issues and pushes to `main`. Unlike atp-siege,
pushing to `main` does **not** auto-deploy anything — Oligolia only builds installers on a
deliberate version-tag push. The remote agent should push fixes to `main` and leave
release-cutting to a deliberate, separate step (see `CLAUDE.md` rule 3).

**Source of truth: GitHub issues** (`gh issue list -R moonsoup/oligolia`).

---

## Remote server quick reference

Run from `~/Software/bio_battle/` — the control scripts are shared across projects:

| What | Command |
|---|---|
| See what it's doing | `bash scripts/remote_status.sh` |
| Send it a message | `bash scripts/remote_send.sh "your message"` |
| Health check (stuck?) | `python3 scripts/loop_health_check.py --tmux` |
| Switch remote to another project | `bash scripts/remote_switch.sh <project>` |

---

## Step 1 — ORIENT

```bash
gh issue list -R moonsoup/oligolia --state open --limit 30
```

Look for issues updated recently (remote likely just pushed a fix). Confirm with:

```bash
gh issue view <N> -R moonsoup/oligolia
```

---

## Step 2 — VERIFY RECENT FIXES

There is no live site to re-test against. Verification means:

- **Backend fix:** `cd backend && .venv/bin/python -m pytest tests/ -q` — all tests pass,
  and ideally a new/updated test specifically covers the fixed behavior.
- **GUI fix:** drive the actual widget headlessly and screenshot it —
  `QT_QPA_PLATFORM=offscreen`, build the widget in a `QApplication`, exercise it via its Qt
  API, `widget.grab().save(...)`, assert on the resulting state. See `gui/panels/` for the
  widget under test; there's no existing automated harness for this beyond writing a
  one-off verification script per fix (don't commit these — they're throwaway checks).
- **Lint:** `ruff check backend/ gui/` must pass for any touched file.

**If verified:**

```bash
gh issue comment <N> -R moonsoup/oligolia --body "## ✅ Verified — $(date +%Y-%m-%d)

Confirmed by: <pytest output / screenshot description>.
Closing."

gh issue close <N> -R moonsoup/oligolia
```

**If NOT resolved:** comment with the failing output and keep the issue open.

---

## Step 3 — DISCOVER NEW ISSUES

```bash
bash .claude/qa/run_all.sh --dry-run   # review findings first
bash .claude/qa/run_all.sh             # files issues (deduped against open issues)
```

This drives the live backend with real NCBI gene data (BRCA1, TP53, HBB, CFTR, EGFR),
diffs expected vs actual output, and traces failures to file:line before filing.

---

## Step 4 — HEALTH CHECK

```bash
python3 scripts/loop_health_check.py --tmux --ping
```

(Run from `~/Software/bio_battle/` — shared script, parameterized by `$CONTAINER`/`$SESSION`/`$VPS`.)

- Exit 0 → remote progressing
- Exit 1 → same issue >10 min; `--ping` auto-sends a kick message
- Exit 2 → VPS unreachable
- Exit 3 → remote is logged out (OAuth token expired/invalid) — see below

**If tmux shows an idle `❯` prompt:**

```bash
bash scripts/remote_send.sh "Resume working the issue queue. gh issue list -R moonsoup/oligolia --state open. Push to main after each fix — do not bump version.py or tag a release."
```

**If exit code is 3 (logged out):**

Copying the local machine's OAuth credential file into the container does **not**
work — tokens are validated server-side and a straight file copy comes back
"Not logged in" even after a fresh restart. The remote needs its own real
`/login` completed through a browser. Claude Code driving this loop should use
`claude-in-chrome` to do it:

```bash
python3 scripts/loop_health_check.py --check-auth
```

This sends `/login`, selects the subscription method, and prints the OAuth
`claude.com/cai/oauth/authorize?...` URL. Then:

1. Open that URL in Chrome (`mcp__claude-in-chrome__navigate`).
2. **Check the signed-in account shown at the bottom of the authorize page
   before clicking anything.** It must be `moonsoup@gmail.com` — the project
   owner's account, matching the local "`moonsoup@gmail.com's Organization`"
   banner Claude Code already shows here. Granting OAuth is an explicit-permission
   action — if the browser is signed into a different account, stop and tell
   the user rather than authorizing the wrong identity.
3. Click **Authorize** — this redirects to `platform.claude.com/oauth/code/callback`
   and displays an authentication code.
4. Read that code off the page (`mcp__claude-in-chrome__get_page_text` or a
   screenshot) and complete sign-in:
   ```bash
   python3 scripts/loop_health_check.py --complete-auth "<code>"
   ```
5. Re-send the resume message (same command as the idle-prompt case above).

Never print the auth code or credential file contents in terminal output —
pipe/relay them directly between tools instead of echoing them.

---

## Rules

- **Never close without running the matching verification first** — only test/screenshot evidence counts.
- **Never bump `version.py` or push a tag** from the remote loop — releases are a deliberate, separate step the human takes.
- **Backend-only issues with passing tests** may close directly; GUI issues need the
  screenshot-verification step above before closing.
