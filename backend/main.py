"""Oligolia — Gene Editing & Viewing Platform — FastAPI backend."""

from contextlib import asynccontextmanager
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fastapi import FastAPI

from version import VERSION
from fastapi.middleware.cors import CORSMiddleware

from .routers import (
    sequences_router, databases_router, files_router,
    alignment_router, crispr_router, variants_router,
    primers_router, pathways_router, analysis_router,
    cloning_router, structure_router,
)


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup
    yield
    # Shutdown


app = FastAPI(
    title="Oligolia",
    description=(
        "Advanced gene editing and viewing platform. "
        "Integrates NCBI, Ensembl, UniProt, KEGG, Reactome, gnomAD, STRING, and PDB. "
        "Supports FASTA, FASTQ, GenBank, EMBL, GFF3, GTF, VCF, BED formats. "
        "Features: sequence editing, CRISPR design, MSA, variant annotation, pathway analysis."
    ),
    version=VERSION,
    lifespan=lifespan,
)

# Permissive CORS is for local development and the QA pipeline only. The shipped
# app makes no HTTP calls to this server at all — the GUI imports the router
# functions and calls them in-process, so nothing here is reachable from a
# browser unless someone runs run_backend.py deliberately (#71).
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(sequences_router)
app.include_router(databases_router)
app.include_router(files_router)
app.include_router(alignment_router)
app.include_router(crispr_router)
app.include_router(variants_router)
app.include_router(primers_router)
app.include_router(pathways_router)
app.include_router(analysis_router)
app.include_router(cloning_router)
app.include_router(structure_router)


@app.get("/")
def root() -> dict:
    return {
        "name": "Oligolia Gene Platform",
        "version": "0.1.0",
        "endpoints": [
            "/sequences", "/databases", "/files",
            "/alignment", "/crispr", "/variants",
            "/primers", "/pathways", "/structure",
            "/docs", "/redoc",
        ],
    }


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}
