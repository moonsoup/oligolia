# Oligolia — Architecture

A free, offline desktop bioinformatics app. One Python process: a PyQt6 GUI that
imports FastAPI router functions and calls them **in-process**.

> Rewritten 2026-09-18. The previous version described a Tauri + React + Vite app
> with a `src-tauri/` Rust shell talking to a FastAPI server over HTTP. None of
> that exists any more, and a document describing an architecture the code does
> not have is worse than no document — it misleads every reader, human or agent
> (#71).

## There is no server

This is the part that changes how everything else reads.

`oligolia.py` never starts uvicorn. The GUI panels import the router functions
directly and call them on a worker thread:

```python
# gui/panels/crispr_panel.py
from backend.routers.crispr import design_guides
...
self._worker = Worker(design_guides, req)
```

No socket, no port, no server process, no CORS, no HTTP in the shipped app. The
FastAPI decorators are still there and the Pydantic models still do the
validation, but nothing serves them to a client.

`run_backend.py` *does* start uvicorn on 127.0.0.1:8765 with `reload=True`. That
is for development and for the QA pipeline under `.claude/qa/`, which drives the
app over HTTP on purpose. Users never run it.

## Layout

```
oligolia/
├── oligolia.py          entry point: QApplication, crash guard, MainWindow
├── version.py           single source of VERSION
├── gui/
│   ├── main_window.py   tabs, menus, the update check
│   ├── workers.py       Worker (QThread) + WorkerSlot/worker_busy guards
│   ├── crash_guard.py   sys.excepthook so a slot exception is a dialog, not exit -6
│   ├── updater.py       GitHub Releases check, patch download + atomic apply
│   └── panels/          one module per tab; each calls backend routers in-process
├── backend/
│   ├── main.py          the FastAPI app — used by run_backend.py and the tests
│   ├── routers/         the endpoints, and where the logic lives
│   ├── services/        external API clients (NCBI, Ensembl, UniProt, …) + helpers
│   ├── models/          Pydantic request/response models
│   ├── formats/         FASTA, FASTQ, GenBank, EMBL, VCF, GFF3, SnapGene, vendor export
│   ├── workflow/         the step engine behind the Workflow tab
│   └── tests/           427 tests (3 skipped without an external aligner)
├── structure_viewer/    OPTIONAL separate app — 3D viewer, its own release cadence
├── frontend/            DEAD. See frontend/DEPRECATED.md
├── .claude/qa/          QA pipeline: scout → runner → analyst → reporter
└── docs/                the marketing site (GitHub Pages) and these notes
```

## Stack

| Layer | Technology |
|---|---|
| GUI | PyQt6 |
| Logic | FastAPI routers + Pydantic models, called in-process |
| Bio | Biopython 1.85 — also the reference implementation the oracle tests check against |
| Packaging | PyInstaller → DMG (macOS, arm64) / Inno Setup (Windows) / AppImage (Linux) |
| CI | `.github/workflows/ci.yml` — lint + backend + offscreen-Qt GUI, gating releases |

## Threading

Anything slow runs on a `Worker` (a `QThread` wrapping a callable). Panels guard
against starting a second job while one is live — rebinding the attribute used to
drop the last reference to a running QThread, which Qt aborts the process over.
Use `WorkerSlot` in new code; `worker_busy(self, "_worker")` is the one-line guard
the existing panels use.

## What is GUI-reachable

Not everything in `backend/` has a tab. ClinVar/gnomAD variant annotation, the
protein property calculator and the ORF finder are backend-only today, and
Synthesis Export is present but gated. Worth knowing before assuming a feature is
user-visible because an endpoint exists.

## Offline vs network

The bioinformatics is local: alignment, CRISPR design, primer design, restriction
mapping, format conversion, ORFs, repeats, protein properties. The network is
used only for database search, variant annotation, structure lookup and the
update check — and 5 of the 12 databases declared in `Database` are actually wired
into `/databases/search`.

Multiple sequence alignment needs an external aligner (MUSCLE or ClustalW) on
PATH. Neither is bundled, so MSA returns a 503 that says so rather than an
approximation — see #76 for the replacement.
