# Toolchain scan — 2026-09-17

What the SPIndle toolchain says about Oligolia, what it got right, and what it cannot
currently be trusted on. Recorded so the next session starts from measurements rather than
from a re-scan.

**Toolchain:** SPIndlebox 1.6.0 (with the plugin contract), Ziggurat 0.2.0 as a SPIndlebox
plugin, stop-guessing 0.5.4. All three were repaired in the same session; the repairs are why
some of these numbers are new.

## Architecture: clean by these checks

`ziggurat report . --full` — **126 source files read, no structural findings.** Three
change-coupling findings, all adjudicated as benign:

| pair | ratio | verdict |
|---|---|---|
| `README.md` ↔ `docs/index.html` | 86% (6 of 7/8 commits) | the site restates the README; expected |
| `docs/css/style.css` ↔ `docs/index.html` | 100% (4 commits) | a page and its stylesheet; expected |
| `gui/update_dialog.py` ↔ `gui/updater.py` | 100% (4 commits) | a dialog and the updater it drives; they genuinely belong together |

No entry-point sprawl, no dynamic loading, no scattered constants or paths, no
sibling-from-global, no singleton bottleneck. Worth stating plainly because the checker was
materially blind a week ago: it could not see `os.makedirs`-style path calls, skipped
`build/`, and matched folder names against the absolute path. Those are fixed
(moonsoup/ziggurat#12, #14, #15, #17, #26, #29), which is why this run reads 126 files rather
than 122.

## Index: fresh

`spindlebox index` → **2494 items, 225 signature classes, 2550 ctx keys**; languages Python,
JavaScript, TypeScript. `spindlebox stale` → up to date, 136 files verified.

## Typing: one real hotspot, and a measurement to distrust

`spindlebox report typing-health` reports **79.9% untyped over 5931 type slots** — but
**3720 of those slots (63%) are `structure_viewer/assets/3Dmol-min.js`**, a vendored minified
third-party bundle. Excluding it, Oligolia is roughly **46% untyped**, and the real hotspots
are its own code:

- `gui/panels/sequence_panel.py` — 64 of 145 slots untyped (the largest panel, and the one
  with the most open correctness issues: #59, #67)
- `backend/tests/test_updater.py` — 64 of 64

Reported to the toolchain as moonsoup/spindlebox#33; until it lands, read the headline
percentage as wrong and the per-file numbers as right.

## Gap analysis: not usable on this project yet

`spindlebox gaps` produced **3005 findings, none actionable**: FastAPI route handlers and Qt
slots counted as dead code (1076), Qt signal handlers as missing ctx keys (905), imported
stdlib names like `StringIO` and methods like `upper` as unresolvable calls (793), and every
no-argument function grouped as a near-duplicate (231). Filed as
moonsoup/spindlebox#35. **Do not spend time triaging that output** until it is fixed.

## Custody: verified

`stop-guessing ledger verify --path .stop-guessing/ledger/custody.jsonl --keyfile <key>` →
**PASS: 1033 records, chain intact and verified under its key.** This was impossible before
today (`--path` recursed until the interpreter gave up; mEllergrace/stop-guessing#94).

## Test and lint state: restored

The documented commands work again. `make venv` had never been run on this machine, so
`backend/.venv` did not exist and both the CLAUDE.md command and the Makefile's `PYTHON`
pointed at a missing interpreter (#70).

```
cd backend && .venv/bin/python -m pytest tests/ -q     213 passed
QT_QPA_PLATFORM=offscreen backend/.venv/bin/python -m pytest \
    gui/test_history.py gui/panels/test_feature_colors.py -q     11 passed
ruff check backend/ gui/                               All checks passed
```

Note for a fresh clone: **run `make venv` first.** Nothing in the repo does it for you, and
no CI job runs the suite at all (#70).

## What to work on next

The scan found no new architectural problems, so the backlog is the correctness work already
filed. In priority order:

1. **#53 — every SpCas9 guide is shifted one base and includes a PAM base.** Confirmed by
   running the code; wrong guides reach the results table and the TSV export, and the existing
   test bakes the bug in. This is the one that costs users real experiments.
2. **#56 — pairwise alignment materialises every co-optimal alignment** (MemoryError on two
   unrelated 500-nt sequences), and it runs on the UI thread.
3. **#54, #55 — the app aborts** on a second concurrent job, and on any unhandled slot
   exception.
4. **#48 — the updater applies unverified code** over the live app, with a bypassable
   traversal guard.
5. **#69, #70 — the QA pipeline cannot run** (gitignored, wrong port) and no CI job runs the
   tests, so nothing above would be caught automatically.

#53 first: it is a one-character regex fix plus its test, and it is currently shipping wrong
answers to anyone who uses the CRISPR panel.
