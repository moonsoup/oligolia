"""Pathway and protein interaction endpoints."""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from ..services import ReactomeClient, KEGGClient, STRINGClient

router = APIRouter(prefix="/pathways", tags=["pathways"])

_reactome = ReactomeClient()
_kegg = KEGGClient()
_string = STRINGClient()


@router.get("/reactome/search")
def reactome_search(query: str, species: str = "Homo sapiens", page_size: int = 20) -> dict:
    try:
        return _reactome.search_query(query, species=species, page_size=page_size)
    except Exception as e:
        raise HTTPException(502, f"Reactome search failed: {e}")


@router.get("/reactome/gene/{gene_id}")
def reactome_gene_pathways(gene_id: str, species: str = "Homo sapiens") -> list[dict]:
    """Get all Reactome pathways containing a given gene (Ensembl/UniProt ID)."""
    try:
        return _reactome.pathways_for_gene(gene_id, species=species)
    except Exception as e:
        raise HTTPException(502, f"Reactome pathway lookup failed: {e}")


@router.get("/reactome/enrichment")
def reactome_enrichment(identifiers: str, species: str = "Homo sapiens") -> dict:
    """Pathway enrichment analysis. Pass comma-separated gene/protein identifiers."""
    ids = [i.strip() for i in identifiers.split(",") if i.strip()]
    if not ids:
        raise HTTPException(400, "No identifiers provided")
    try:
        return _reactome.mapping_analysis(ids, species=species)
    except Exception as e:
        raise HTTPException(502, f"Reactome analysis failed: {e}")


@router.get("/kegg/search")
def kegg_search(query: str, database: str = "pathway") -> list[dict]:
    try:
        results = _kegg.search(database, query)
        return [{"id": r[0], "description": r[1]} for r in results]
    except Exception as e:
        raise HTTPException(502, f"KEGG search failed: {e}")


@router.get("/kegg/pathway/{pathway_id}/genes")
def kegg_pathway_genes(pathway_id: str) -> dict:
    try:
        genes = _kegg.genes_in_pathway(pathway_id)
        return {"pathway_id": pathway_id, "genes": genes, "count": len(genes)}
    except Exception as e:
        raise HTTPException(502, f"KEGG pathway genes failed: {e}")


class StringNetworkRequest(BaseModel):
    proteins: list[str]
    species: int = 9606
    required_score: int = 400


@router.post("/string/network")
def string_network(req: StringNetworkRequest) -> list[dict]:
    """Get STRING protein-protein interaction network for a list of proteins."""
    try:
        return _string.network(req.proteins, species=req.species, required_score=req.required_score)
    except Exception as e:
        raise HTTPException(502, f"STRING network failed: {e}")


@router.get("/string/partners/{protein}")
def string_partners(protein: str, species: int = 9606, limit: int = 10, min_score: int = 700) -> list[dict]:
    try:
        return _string.interaction_partners(protein, species=species, limit=limit, required_score=min_score)
    except Exception as e:
        raise HTTPException(502, f"STRING partners failed: {e}")


@router.get("/string/enrichment")
def string_enrichment(proteins: str, species: int = 9606) -> list[dict]:
    ids = [p.strip() for p in proteins.split(",") if p.strip()]
    try:
        return _string.enrichment(ids, species=species)
    except Exception as e:
        raise HTTPException(502, f"STRING enrichment failed: {e}")
