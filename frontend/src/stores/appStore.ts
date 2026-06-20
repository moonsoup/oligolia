import { create } from "zustand";
import type { Sequence, Variant } from "../types";

type Panel = "sequences" | "search" | "editor" | "crispr" | "variants" | "alignment" | "primers" | "pathways" | "genome";

interface AppStore {
  // Active panel
  activePanel: Panel;
  setActivePanel: (p: Panel) => void;

  // Active sequence (for editor)
  activeSequenceId: string | null;
  setActiveSequenceId: (id: string | null) => void;

  // Loaded sequences
  sequences: Sequence[];
  setSequences: (seqs: Sequence[]) => void;
  addSequence: (seq: Sequence) => void;
  removeSequence: (id: string) => void;

  // Variants
  variants: Variant[];
  setVariants: (v: Variant[]) => void;

  // File upload result
  lastUploadedSequences: Sequence[];
  setLastUploadedSequences: (seqs: Sequence[]) => void;

  // Backend URL (for Tauri sidecar)
  backendUrl: string;
  setBackendUrl: (url: string) => void;
}

export const useAppStore = create<AppStore>((set) => ({
  activePanel: "sequences",
  setActivePanel: (p) => set({ activePanel: p }),

  activeSequenceId: null,
  setActiveSequenceId: (id) => set({ activeSequenceId: id }),

  sequences: [],
  setSequences: (sequences) => set({ sequences }),
  addSequence: (seq) => set((s) => ({ sequences: [...s.sequences, seq] })),
  removeSequence: (id) =>
    set((s) => ({ sequences: s.sequences.filter((seq) => seq.id !== id) })),

  variants: [],
  setVariants: (variants) => set({ variants }),

  lastUploadedSequences: [],
  setLastUploadedSequences: (seqs) => set({ lastUploadedSequences: seqs }),

  backendUrl: "http://127.0.0.1:8765",
  setBackendUrl: (backendUrl) => set({ backendUrl }),
}));
