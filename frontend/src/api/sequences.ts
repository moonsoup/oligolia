import { api } from "./client";
import type { Sequence, SequenceEditRequest, SequenceEditResult } from "../types";

export const sequences = {
  list: () => api.get<Sequence[]>("/sequences/").then((r) => r.data),
  get: (id: string) => api.get<Sequence>(`/sequences/${id}`).then((r) => r.data),
  add: (seq: Sequence) => api.post<Sequence>("/sequences/", seq).then((r) => r.data),
  delete: (id: string) => api.delete(`/sequences/${id}`),

  edit: (id: string, req: SequenceEditRequest) =>
    api.post<SequenceEditResult>(`/sequences/${id}/edit`, req).then((r) => r.data),

  gcContent: (id: string) =>
    api.get<{ seq_id: string; length: number; gc_content: number }>(`/sequences/${id}/gc_content`).then((r) => r.data),

  codonUsage: (id: string) =>
    api.get(`/sequences/${id}/codon_usage`).then((r) => r.data),

  findMotif: (id: string, motif: string) =>
    api.get(`/sequences/${id}/find_motif`, { params: { motif } }).then((r) => r.data),
};
