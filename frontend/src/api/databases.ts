import { api } from "./client";
import type { SearchRequest, SearchResponse } from "../types";

export const databases = {
  search: (req: SearchRequest) =>
    api.post<SearchResponse>("/databases/search", req).then((r) => r.data),

  fetchNCBISequence: (geneId: string, format: "fasta" | "gb" = "fasta") =>
    api.get<{ content: string }>(`/databases/ncbi/gene/${geneId}/sequence`, { params: { format } }).then((r) => r.data),

  fetchEnsemblSequence: (ensemblId: string, seqType = "cdna") =>
    api.get(`/databases/ensembl/${ensemblId}/sequence`, { params: { seq_type: seqType } }).then((r) => r.data),

  clinvarSearch: (gene: string, maxResults = 20) =>
    api.get(`/databases/ncbi/clinvar/${gene}`, { params: { max_results: maxResults } }).then((r) => r.data),

  blast: (sequence: string, program = "blastn", database = "nt", maxHits = 10) =>
    api.post("/databases/blast", null, {
      params: { sequence, program, database, max_hits: maxHits },
    }).then((r) => r.data),

  keggPathway: (pathwayId: string) =>
    api.get(`/databases/kegg/pathway/${pathwayId}`).then((r) => r.data),
};
