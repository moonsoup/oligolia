"""Multi-database search and fetch endpoints."""

from fastapi import APIRouter, HTTPException, Query
from ..models.search import SearchRequest, SearchResponse, SearchResult, Database
from ..services import NCBIClient, EnsemblClient, UniProtClient, KEGGClient

router = APIRouter(prefix="/databases", tags=["databases"])

_ncbi = NCBIClient()
_ensembl = EnsemblClient()
_uniprot = UniProtClient()
_kegg = KEGGClient()


@router.post("/search", response_model=SearchResponse)
def search(req: SearchRequest) -> SearchResponse:
    results: list[SearchResult] = []
    errors: dict[str, str] = {}

    for db in req.databases:
        try:
            if db == Database.NCBI_GENE:
                records = _ncbi.search_genes(req.query, organism=req.species, max_results=req.max_results)
                for r in records:
                    results.append(SearchResult(
                        id=str(r.get("uid", "")),
                        name=r.get("name", ""),
                        description=r.get("description", ""),
                        database=db,
                        organism=r.get("organism", {}).get("scientificname", req.species),
                        accession=str(r.get("uid", "")),
                        url=f"https://www.ncbi.nlm.nih.gov/gene/{r.get('uid', '')}",
                    ))

            elif db == Database.NCBI_NUCLEOTIDE:
                result = _ncbi.esearch("nucleotide", req.query, retmax=req.max_results)
                for uid in result.get("esearchresult", {}).get("idlist", []):
                    results.append(SearchResult(
                        id=uid, name=uid, description=f"Nucleotide record {uid}",
                        database=db, accession=uid,
                        url=f"https://www.ncbi.nlm.nih.gov/nuccore/{uid}",
                    ))

            elif db == Database.ENSEMBL:
                try:
                    gene = _ensembl.lookup_symbol(req.species.lower().replace(" ", "_"), req.query)
                    if gene:
                        results.append(SearchResult(
                            id=gene.get("id", ""),
                            name=gene.get("display_name", req.query),
                            description=gene.get("description", ""),
                            database=db,
                            organism=gene.get("species", req.species),
                            accession=gene.get("id", ""),
                            url=f"https://www.ensembl.org/id/{gene.get('id', '')}",
                            extra={
                                "biotype": gene.get("biotype", ""),
                                "assembly": gene.get("assembly_name", ""),
                                "location": f"{gene.get('seq_region_name', '')}:{gene.get('start', '')}-{gene.get('end', '')}",
                            },
                        ))
                except Exception:
                    pass  # symbol not found is expected

            elif db == Database.UNIPROT:
                data = _uniprot.search_by_gene(req.query, organism=req.species, size=req.max_results)
                for entry in data.get("results", []):
                    acc = entry.get("primaryAccession", "")
                    gene_names = entry.get("genes", [{}])
                    gene_name = gene_names[0].get("geneName", {}).get("value", req.query) if gene_names else req.query
                    results.append(SearchResult(
                        id=acc, name=gene_name,
                        description=entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
                        database=db, organism=req.species, accession=acc,
                        url=f"https://www.uniprot.org/uniprotkb/{acc}",
                    ))

            elif db == Database.KEGG:
                matches = _kegg.search("genes", req.query)
                for kegg_id, desc in matches[:req.max_results]:
                    results.append(SearchResult(
                        id=kegg_id, name=kegg_id, description=desc,
                        database=db, accession=kegg_id,
                        url=f"https://www.kegg.jp/entry/{kegg_id}",
                    ))

        except Exception as e:
            errors[db.value] = str(e)

    return SearchResponse(
        query=req.query,
        total=len(results),
        results=results,
        databases_searched=req.databases,
        errors=errors,
    )


@router.get("/ncbi/gene/{gene_id}/sequence")
def fetch_ncbi_gene_sequence(gene_id: str, format: str = Query("fasta", enum=["fasta", "gb"])) -> dict:
    """Fetch mRNA sequence for an NCBI Gene UID (converts Gene UID → RefSeq nuccore first)."""
    try:
        text = _ncbi.fetch_fasta_for_gene(gene_id)
        if format == "gb":
            # Re-fetch in GenBank format using the resolved nuccore ID
            nuccore_ids = _ncbi.gene_uid_to_refseq_rna(gene_id)
            if nuccore_ids:
                text = _ncbi.fetch_nucleotide(nuccore_ids[0])
        return {"gene_id": gene_id, "format": format, "content": text}
    except Exception as e:
        raise HTTPException(502, f"NCBI fetch failed: {e}")


@router.get("/ensembl/{ensembl_id}/sequence")
def fetch_ensembl_sequence(ensembl_id: str, seq_type: str = Query("cdna", enum=["genomic", "cdna", "cds", "protein"])) -> dict:
    try:
        data = _ensembl.sequence_id(ensembl_id, seq_type=seq_type)
        return {"ensembl_id": ensembl_id, "seq_type": seq_type, "sequence": data.get("seq", ""), "desc": data.get("desc", "")}
    except Exception as e:
        raise HTTPException(502, f"Ensembl fetch failed: {e}")


@router.get("/ncbi/clinvar/{gene}")
def clinvar_search(gene: str, max_results: int = Query(20, ge=1, le=100)) -> dict:
    try:
        records = _ncbi.search_clinvar(gene, max_results)
        return {"gene": gene, "total": len(records), "records": records}
    except Exception as e:
        raise HTTPException(502, f"ClinVar search failed: {e}")


@router.post("/blast")
def blast_search(sequence: str, program: str = "blastn", database: str = "nt", max_hits: int = 10) -> dict:
    try:
        return _ncbi.blast_search(sequence, program=program, database=database, max_hits=max_hits)
    except TimeoutError:
        raise HTTPException(504, "BLAST search timed out")
    except Exception as e:
        raise HTTPException(502, f"BLAST failed: {e}")


@router.get("/kegg/pathway/{pathway_id}")
def kegg_pathway(pathway_id: str) -> dict:
    try:
        text = _kegg.get_pathway(pathway_id)
        return {"pathway_id": pathway_id, "content": text}
    except Exception as e:
        raise HTTPException(502, f"KEGG fetch failed: {e}")
