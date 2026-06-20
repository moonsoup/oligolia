// ── Sequences ────────────────────────────────────────────────────────────────

export type MoleculeType = "DNA" | "RNA" | "PROTEIN" | "UNKNOWN";

export interface Annotation {
  feature_type: string;
  start: number;
  end: number;
  strand: "+" | "-" | ".";
  qualifiers: Record<string, unknown>;
}

export interface Sequence {
  id: string;
  name: string;
  description: string;
  seq: string;
  molecule_type: MoleculeType;
  annotations: Annotation[];
  source_db: string;
  accession: string;
  length: number;
}

export interface SequenceEditRequest {
  operation:
    | "insert"
    | "delete"
    | "replace"
    | "reverse_complement"
    | "complement"
    | "translate"
    | "transcribe"
    | "back_transcribe";
  position?: number;
  end_position?: number;
  insert_seq?: string;
  replacement?: string;
}

export interface SequenceEditResult {
  original_id: string;
  operation: string;
  result_seq: string;
  diff_start?: number;
  diff_end?: number;
  message: string;
}

// ── Search ───────────────────────────────────────────────────────────────────

export type Database =
  | "ncbi_gene"
  | "ncbi_nucleotide"
  | "ncbi_protein"
  | "ncbi_snp"
  | "ncbi_clinvar"
  | "ensembl"
  | "uniprot"
  | "kegg"
  | "reactome"
  | "gnomad"
  | "string"
  | "pdb";

export interface SearchRequest {
  query: string;
  databases: Database[];
  species: string;
  max_results: number;
}

export interface SearchResult {
  id: string;
  name: string;
  description: string;
  database: Database;
  organism: string;
  accession: string;
  url: string;
  extra: Record<string, string>;
}

export interface SearchResponse {
  query: string;
  total: number;
  results: SearchResult[];
  databases_searched: Database[];
  errors: Record<string, string>;
}

// ── Variants ─────────────────────────────────────────────────────────────────

export type VariantType = "SNP" | "INDEL" | "DEL" | "INS" | "CNV" | "SV" | "UNKNOWN";

export interface Variant {
  chrom: string;
  pos: number;
  ref: string;
  alt: string[];
  id: string;
  qual?: number;
  filter: string[];
  info: Record<string, unknown>;
  variant_type: VariantType;
  gene: string;
  clinical_significance?: string;
  allele_frequency?: number;
  gnomad_af?: number;
}

// ── CRISPR ───────────────────────────────────────────────────────────────────

export type CasType = "SpCas9" | "SpCas9-HF1" | "AsCas12a" | "LwaCas13a";

export interface GuideRNA {
  sequence: string;
  pam: string;
  position: number;
  strand: "+" | "-";
  gc_content: number;
  on_target_score?: number;
  off_target_count?: number;
  efficiency_score?: number;
}

export interface CRISPRDesignRequest {
  target_sequence: string;
  cas_type: CasType;
  guide_length: number;
  max_guides: number;
  check_off_targets: boolean;
}

export interface CRISPRDesignResponse {
  target_length: number;
  cas_type: CasType;
  guides: GuideRNA[];
  total_candidates: number;
}

// ── Alignment ─────────────────────────────────────────────────────────────────

export interface AlignmentResult {
  score: number;
  aligned_seq1: string;
  aligned_seq2: string;
  identity: number;
  similarity: number;
  gaps: number;
  alignment_length: number;
}

// ── Primers ───────────────────────────────────────────────────────────────────

export interface Primer {
  sequence: string;
  position: number;
  length: number;
  tm: number;
  gc_content: number;
  direction: "forward" | "reverse";
}

export interface PrimerPair {
  forward: Primer;
  reverse: Primer;
  product_size: number;
  penalty: number;
}

export interface RestrictionSite {
  enzyme: string;
  cut_pattern: string;
  positions: number[];
  count: number;
}
