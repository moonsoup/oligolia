# Overnight run — 2026-09-17 into 2026-09-18

What was done, what it cost, and what is left. The GitHub issues carry the
evidence; this is the map.

**25 commits on `main` past `v0.5.2`. 24 issues closed with verified evidence, 8
new ones filed. 432 backend + 88 GUI tests, up from 213 + 11.**

## The gate that did not exist

`.github/workflows/ci.yml` runs lint, the backend suite and the offscreen-Qt GUI
suite on every push and PR, and `release.yml` will not build anything until it
passes (#70). Before this, nothing ran `pytest` or `ruff` in CI and releases were
built from unverified code. Its first run on GitHub matched local exactly.

`make check` runs the same three steps locally. `make test` is still the backend
suite alone.

## Closed, with the defect reproduced first in every case

| # | What was wrong |
|---|---|
| 48 | the updater could apply a truncated patch over the live bundle and brick it; traversal guard was a string prefix, so `Oligolia.app.evil` passed |
| 50 | primer Tm was the 1989 Wallace rule under a nearest-neighbour name — 97 of 120 primers >3 °C out |
| 52 | no 3′ primer-dimer check between the pair, which no per-primer filter can catch |
| 54 | a second concurrent job aborted the process (exit 134); nine call sites |
| 55 | no `sys.excepthook`, so any unhandled slot exception killed the app |
| 56 | pairwise alignment materialised every co-optimal alignment — OverflowError at 500 nt |
| 57 | GenBank round-trip collapsed `join()` to outer bounds and stringified qualifier lists |
| 58 | MSA right-padded the input and called it an alignment |
| 59 | edits never shifted annotations; `Replace` with End=0 grew the sequence |
| 60 | Optimize Codons mutated the sequence then called a method that does not exist |
| 61 | ClinVar applied one record's significance to every variant in the gene, matched by substring |
| 62 | hairpin filter rejected 56% of candidates; reverse primers skipped the forward rules |
| 63 | four protein properties disagreed with ProtParam, all in a misleading direction |
| 64 | Cas12a scanned one strand, Cas13 dropped its last window, `guide_length` ignored |
| 65 | one microsatellite came back as 50 fragments; a real repeat was dropped |
| 66 | SASA computed across all NMR models, so they occluded each other |
| 67 | `<DEL>` typed as an insertion; `.embl` files loaded as nothing |
| 68 | BLAST asked for JSON2, which NCBI returns as a zip — every search timed out |
| 69 | the QA pipeline was gitignored, on the wrong port, with a stale endpoint map |
| 70 | two workflows raced on every tag; no CI ran the tests |
| 73 | QA scout substituted HBB for any gene whose fetch failed, and the run passed |
| 74 | QA reporter bypassed `gh_issue.py` and deduped on a two-word title overlap |
| 77 | the alignment view's match bars wrapped out of step with the sequences |
| 80 | the QA pipeline's inter-phase state wrote to a file that does not exist |

## Three things worth knowing beyond the list

**The QA pipeline runs, and its first real run would have filed 28 false issues.**
Every one was a defect in its own corpus — a `cas_type` that is not in the enum, a
missing trailing slash, body params that are query params, field names that never
existed. Corrected across four runs: 28 → 7 → 1 → 0. The app passes its own
pipeline. `--dry-run` first, always.

**Codex reviewed the batch read-only and made four hits on my own work.** The
sharpest: my "differential oracle" tests wrapped `Tm_NN` and then compared against
`Tm_NN` — a function to itself. Both oracle files now assert against values
recorded from Biopython 1.85, with a separate test that fails if Biopython's own
tables drift. That mattered beyond the two files: #75 proposes a whole oracle
layer, and "wrap the oracle, then compare to the oracle" would have made it
vacuous by construction. It also asked whether the crash guard had traded a crash
for a dialog loop. It had — 55 dialogs in 400 ms, measured.

**Looking at the app found a bug no assertion would have.** Driving the real
window offscreen and reading the screenshots showed the alignment view's match
bars sitting under the wrong bases (#77): score, identity and gap count all
correct, picture wrong. The layout findings are #78.

## Still open

- **#28, #29** — legal, waiting on attorney review. Untouched.
- **#51** — BLAST specificity in primer design. Unblocked by #68; not built.
- **#71** — the doc corrections are done; the `structview-v*` release, the issue
  digest failing on `KeyError: 'url'`, the placeholder worker URL in
  `docs/js/integrate.js`, and rule 5 contradicting `release.yml` are not.
- **#72** — the stray `main` tag. Confirmed safe to delete; two commands left for
  the owner, because deleting a remote tag is not reversible from here.
- **#75** — the oracle layer. Tm and protein properties landed; the coordinate and
  guide-re-derivation checks exist as regular tests rather than a layer.
- **#76** — MSA has no working default now that the fake one is gone.
- **#78** — layout: the 10th tab is clipped at 950 px height.
- **#79** — Tm is pinned at `Mg=0` and 25 nM primer, which is thermodynamically
  correct for a buffer nobody runs. Needs a scientific decision, not a patch.

## A version bump is owed

25 commits past `v0.5.2`, including protein structure prediction. Per CLAUDE.md
rule 3 that decision is the interactive session's, not mine to take unprompted —
and #79 and #76 are both open questions about behaviour users will see. The CI
gate means a release build will not start until lint and both suites pass.

## One process failure worth recording

I pushed `1152087` after running only the new test file, and two tests were red on
`main` for one commit. Fixed in `b08fcda`. Running the full suite before every
commit is not optional and I skipped it.
