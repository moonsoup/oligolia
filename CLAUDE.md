# Oligolia — CLAUDE.md

**Tagline:** Gene editing for everyone — free, offline, desktop.
**Repo:** https://github.com/moonsoup/oligolia (public)

## What This Is

A free, open-source desktop bioinformatics platform: PyQt6 GUI + an embedded FastAPI
backend (not a hosted service). Search NCBI/Ensembl/UniProt/etc., design CRISPR guides,
align sequences, annotate variants — all offline, packaged as a native installer.

## Stack

| Layer | Tech |
|---|---|
| GUI | PyQt6 (`gui/`) |
| Backend | FastAPI, embedded and launched by the GUI process (`backend/`) |
| Bio | Biopython |
| Packaging | PyInstaller → DMG (macOS) / Inno Setup (Windows) / AppImage (Linux) |
| Distribution | GitHub Releases — **built on version-tag push only**, not on every push to main |
| Repo | moonsoup/oligolia (public) |

## Running / Testing Locally

```bash
python oligolia.py                                  # launch the GUI
cd backend && .venv/bin/python -m pytest tests/ -q   # backend test suite
ruff check backend/ gui/                             # lint
```

For GUI changes, verify by driving the real widget headlessly rather than asserting on
internals only: `QT_QPA_PLATFORM=offscreen`, build the widget in a `QApplication`, drive
it via its Qt API (set fields, click buttons), `widget.grab().save(...)` a screenshot, and
assert on the resulting state/text. This is the closest thing this project has to a browser
e2e test, since there is no live website to drive.

## Development Rules

1. **Lint before commit:** `ruff check backend/ gui/` (or the specific changed files). Never commit failing lint.
2. **Push directly to `main`** — solo/small-team project, no PR review ceremony.
3. **Version bumps are batched, not per-fix.** Land fixes on `main` as normal commits with
   no version change. Only bump `version.py`, commit "Bump to X.Y.Z", tag `vX.Y.Z`, and push
   the tag when actually cutting a release — that tag push is what triggers the GitHub
   Actions build/release. Pushing to `main` alone does **not** auto-release or auto-deploy
   (unlike atp-siege's deploy-on-push model).
4. **Release assets need two filenames:** a stable one (`Oligolia-mac.dmg`, used by
   Pages/README, never changes) and a versioned one (`Oligolia-X.Y.Z-mac.dmg`). Every
   release must upload both.
5. After a release tag: run `scripts/make_patch.py` (build patch + manifest), then
   `gh release create vX.Y.Z` with the patch/manifest assets attached.

## QA Pipeline (bug discovery — the equivalent of a live-site playtest agent)

Oligolia is a desktop app, not a website — there's no Playwright/browser specialty here.
`.claude/qa/` is the equivalent: real-NCBI-data-driven backend testing that diffs expected
vs actual output and files GitHub issues automatically.

```bash
bash .claude/qa/run_all.sh --dry-run   # scout → runner → analyst, no issues filed — review first
bash .claude/qa/run_all.sh             # also files issues via reporter.py (deduped against open issues)
```

Pipeline: `scout.py` (builds a corpus from real genes — BRCA1, TP53, HBB, CFTR, EGFR — via
NCBI/Ensembl) → `runner.py` (drives the live backend, `python run_backend.py` must be
running) → `analyst.py` (diffs expected vs actual, traces failures to file:line) →
`reporter.py` (files deduped GitHub issues).

## Synthesis Vendor Export

`backend/formats/synthesis_order.py` generates vendor-formatted order files (GENEWIZ xlsx,
IDT/Twist/Eurofins CSV) for a sequence — never submits a live order. `VENDORS` dict holds
each vendor's column layout and hand-off instructions; `verified=False` means the format
wasn't confirmed against the vendor's own published spec and should be double-checked
before relying on it.

## Remote Dev Loop (Docker Claude on the Hostinger VPS)

A Claude Code instance can run continuously on the VPS inside the shared `claude-remote`
container, working through `gh issue list -R moonsoup/oligolia` autonomously, the same
mechanism used for atp-siege (see `~/Software/bio_battle/CLAUDE.md` for the full
architecture). The container is **shared across projects** via an active-project marker —
see `~/Software/bio_battle/scripts/remote_switch.sh oligolia` to point it here, and
`remote_switch.sh atpsiege` to hand it back.

**Quick reference (run from `~/Software/bio_battle/`, the scripts live there):**

| What | Command |
|---|---|
| See what it's doing | `bash scripts/remote_status.sh` |
| Send it a message | `bash scripts/remote_send.sh "message"` |
| Switch active project | `bash scripts/remote_switch.sh oligolia` |
| Health check (stuck?) | `python3 scripts/loop_health_check.py --tmux` |

**Differences from atp-siege's loop:**
- No live-site deploy-on-push — pushing to `main` just lands the fix. A release (and the
  installers it produces) only happens on a deliberate version-tag push (see rule 3 above).
  The remote agent should push fixes to `main` but **not** cut releases itself.
- DISCOVER step uses `.claude/qa/run_all.sh`, not a Playwright/axe-core playtest agent.
- VERIFY step re-runs `cd backend && .venv/bin/python -m pytest tests/ -q` for backend
  fixes, or the offscreen-Qt screenshot pattern above for GUI fixes — there's no live URL
  to re-test against.

See `.claude/skills/advance-queue/SKILL.md` in this repo for the full loop, adapted from
atp-siege's `/advance-queue`.
