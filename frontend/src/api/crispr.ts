import { api } from "./client";
import type { CRISPRDesignRequest, CRISPRDesignResponse } from "../types";

export const crispr = {
  design: (req: CRISPRDesignRequest) =>
    api.post<CRISPRDesignResponse>("/crispr/design", req).then((r) => r.data),

  scoreGuide: (guideSequence: string) =>
    api.post(`/crispr/score_guide`, null, { params: { guide_sequence: guideSequence } }).then((r) => r.data),
};
