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

Run from `~/Software/Friday_tools/remote_loop/` — the control scripts are shared across projects:

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

(Run from `~/Software/Friday_tools/remote_loop/` — shared script, parameterized by `$CONTAINER`/`$SESSION`/`$VPS`.)

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

## Step 5 — DESIGN REVIEW (every ~4-5 iterations, or after a batch of new UI lands)

Functional verification (Step 2) checks that a fix works, not that the app
still looks and feels professional as new surface area (plasmid map,
feature table, assembly tab, etc.) accumulates. Every few loop iterations —
not every single one — do a design pass:

1. Launch the real app offscreen and screenshot the panels that changed
   recently (`QT_QPA_PLATFORM=offscreen`, real `QApplication`, drive the
   actual widget, `grab().save()` — same pattern as Step 2's GUI
   verification, just looking at the whole panel rather than one fix).
2. Invoke the `frontend-design` skill against those screenshots — its
   guidance on distinctive, intentional visual design and avoiding
   templated-default clutter applies here even though Oligolia is a desktop
   Qt app, not a web frontend: look for crowded/unbalanced layouts,
   inconsistent spacing or type treatment between panels, and controls that
   don't read as purposeful.
3. Prioritize ease of use over feature density — if a panel added this
   cycle makes the app feel more cluttered rather than more capable, that's
   worth a follow-up issue even though nothing is "broken."
4. File findings as normal GitHub issues (`bug` for regressions in
   existing polish, `enhancement` for new cleanup work) rather than fixing
   inline — keeps this step fast and matches the rest of the loop's
   verify-then-file pattern.

---

## Rules

- **Never close without running the matching verification first** — only test/screenshot evidence counts.
- **Never bump `version.py` or push a tag** from the remote loop — releases are a deliberate, separate step the human takes.
- **Backend-only issues with passing tests** may close directly; GUI issues need the
  screenshot-verification step above before closing.
