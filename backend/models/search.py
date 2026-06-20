from enum import Enum
from pydantic import BaseModel, Field


class Database(str, Enum):
    NCBI_GENE = "ncbi_gene"
    NCBI_NUCLEOTIDE = "ncbi_nucleotide"
    NCBI_PROTEIN = "ncbi_protein"
    NCBI_SNP = "ncbi_snp"
    NCBI_CLINVAR = "ncbi_clinvar"
    ENSEMBL = "ensembl"
    UNIPROT = "uniprot"
    KEGG = "kegg"
    REACTOME = "reactome"
    GNOMAD = "gnomad"
    STRING = "string"
    PDB = "pdb"


class SearchRequest(BaseModel):
    query: str
    databases: list[Database] = Field(default_factory=lambda: [Database.NCBI_GENE, Database.ENSEMBL])
    species: str = "human"
    max_results: int = Field(default=20, ge=1, le=200)


class SearchResult(BaseModel):
    id: str
    name: str
    description: str
    database: Database
    organism: str = ""
    accession: str = ""
    url: str = ""
    extra: dict = Field(default_factory=dict)


class SearchResponse(BaseModel):
    query: str
    total: int
    results: list[SearchResult]
    databases_searched: list[Database]
    errors: dict[str, str] = Field(default_factory=dict)
