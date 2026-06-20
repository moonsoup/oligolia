# Oligolia — Architecture

## Overview

Cross-platform gene editing and viewing platform. Backend processes all bioinformatics; frontend provides the UI; Tauri wraps everything into a native desktop app.

```
oligolia/
├── backend/          Python FastAPI — bioinformatics engine + DB API clients
├── frontend/         React + TypeScript + Vite — UI
├── src-tauri/        Rust Tauri shell — packages as native app
├── run_backend.py    Dev entry point for backend
└── docs/
```

## Stack

| Layer | Technology | Why |
|-------|-----------|-----|
| Backend | Python 3.11 + FastAPI + Biopython | Best-in-class bioinformatics ecosystem |
| Frontend | React 18 + TypeScript + Vite | Fast build, strong typing, excellent ecosystem |
| State | Zustand | Simple, no boilerplate |
| Data fetching | TanStack Query | Caching, loading states, retry |
| Desktop shell | Tauri 2.x | 10× smaller than Electron, supports iOS/Android |
| Styling | Tailwind CSS | Rapid dark-mode-first UI |

## Backend Modules

### Services (external API clients)
- `ncbi.py` — NCBI Entrez E-utilities + Datasets v2 + BLAST
- `ensembl.py` — Ensembl REST API v15 (lookup, sequence, VEP, homology)
- `uniprot.py` — UniProt REST API (search, entry, ID mapping)
- `kegg.py` — KEGG REST API (genes, pathways, compounds)
- `reactome.py` — Reactome Content Service (pathways, enrichment)
- `gnomad.py` — gnomAD GraphQL API (allele frequencies, variants)
- `string_db.py` — STRING API v12 (protein interactions, enrichment)
- `pdb.py` — RCSB PDB REST + search (structures, sequences)

### Format Handlers
- `fasta.py` — FASTA + FASTQ read/write
- `genbank.py` — GenBank + EMBL read/write
- `vcf.py` — VCF v4.3 parse/write/stream
- `gff.py` — GFF3 + GTF parse/write

### API Routers
- `/sequences` — CRUD + edit (insert/delete/replace/RC/translate/transcribe)
- `/databases` — multi-DB search, fetch by ID, BLAST
- `/files` — upload/parse/convert/download all formats
- `/alignment` — pairwise (NW/SW) + MSA (MUSCLE/fallback)
- `/crispr` — SpCas9/Cas9-HF/Cas12a/Cas13 guide design + scoring
- `/variants` — VCF annotation with gnomAD/ClinVar
- `/primers` — PCR primer design + restriction enzyme analysis (20 enzymes)
- `/pathways` — Reactome, KEGG, STRING network endpoints

## Databases Integrated

| Database | Access | Data |
|----------|--------|------|
| NCBI Gene | E-utilities REST | Gene records, sequences, summaries |
| NCBI Nucleotide | E-utilities REST | GenBank records, FASTA |
| NCBI ClinVar | E-utilities REST | Clinical variant significance |
| NCBI BLAST | BLAST REST | Sequence similarity search |
| Ensembl | REST API v15 | Gene lookup, VEP, sequences, homology |
| UniProt | REST API | Protein entries, ID mapping |
| KEGG | REST API | Pathways, genes, compounds |
| Reactome | Content Service | Pathway hierarchy, enrichment |
| gnomAD | GraphQL | Allele frequencies, constraint |
| STRING | REST API v12 | Protein-protein interactions |
| RCSB PDB | REST + search | 3D structures |

## File Formats Supported

**Input:** FASTA, FASTQ, GenBank (.gb/.gbk), EMBL, GFF3, GTF, VCF
**Output:** FASTA, FASTQ, GenBank, VCF, GFF3, TSV (guides/primers)

## Cross-Platform Packaging

Development: `python run_backend.py` + `npm run dev` in frontend/

Production (Tauri):
1. `cd backend && pyinstaller --onefile main_sidecar.py`  
2. Place binary in `src-tauri/binaries/`
3. `cd src-tauri && cargo tauri build`

Desktop: Win/Mac/Linux via Tauri 2.x  
Mobile: iOS/Android via Tauri 2.x (future — same React frontend)
