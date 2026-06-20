import { useAppStore } from "../../stores/appStore";
import { sequences as seqApi } from "../../api";
import { useQuery, useQueryClient } from "@tanstack/react-query";

const PANEL_ICONS: Record<string, string> = {
  sequences: "🧬",
  search: "🔍",
  editor: "✏️",
  crispr: "✂️",
  variants: "🔬",
  alignment: "↔️",
  primers: "🔩",
  pathways: "🗺️",
  genome: "🌐",
};

const PANEL_LABELS: Record<string, string> = {
  sequences: "Sequences",
  search: "DB Search",
  editor: "Editor",
  crispr: "CRISPR",
  variants: "Variants",
  alignment: "Alignment",
  primers: "Primers",
  pathways: "Pathways",
  genome: "Genome",
};

export function Sidebar() {
  const { activePanel, setActivePanel, sequences, activeSequenceId, setActiveSequenceId, removeSequence } =
    useAppStore((s) => ({
      activePanel: s.activePanel,
      setActivePanel: s.setActivePanel,
      sequences: s.sequences,
      activeSequenceId: s.activeSequenceId,
      setActiveSequenceId: s.setActiveSequenceId,
      removeSequence: s.removeSequence,
    }));
  const qc = useQueryClient();

  return (
    <div className="w-64 bg-gray-900 border-r border-gray-700 flex flex-col h-full">
      {/* Logo */}
      <div className="p-4 border-b border-gray-700">
        <h1 className="text-xl font-bold text-green-400">Oligolia</h1>
        <p className="text-xs text-gray-500 mt-0.5">Gene Editing Platform</p>
      </div>

      {/* Nav */}
      <nav className="p-2 border-b border-gray-700">
        {Object.keys(PANEL_ICONS).map((panel) => (
          <button
            key={panel}
            onClick={() => setActivePanel(panel as never)}
            className={`w-full flex items-center gap-3 px-3 py-2 rounded-lg text-sm transition-colors mb-0.5 ${
              activePanel === panel
                ? "bg-green-900/40 text-green-300"
                : "text-gray-400 hover:bg-gray-800 hover:text-gray-200"
            }`}
          >
            <span>{PANEL_ICONS[panel]}</span>
            <span>{PANEL_LABELS[panel]}</span>
          </button>
        ))}
      </nav>

      {/* Loaded sequences */}
      <div className="flex-1 overflow-auto p-2">
        <p className="text-xs text-gray-500 px-2 py-1 font-medium uppercase tracking-wide">
          Loaded Sequences
        </p>
        {sequences.length === 0 && (
          <p className="text-xs text-gray-600 px-3 py-2">
            None — search a DB or upload a file.
          </p>
        )}
        {sequences.map((s) => (
          <div
            key={s.id}
            className={`group flex items-center gap-2 px-2 py-1.5 rounded text-xs cursor-pointer transition-colors ${
              activeSequenceId === s.id
                ? "bg-green-900/30 text-green-300"
                : "text-gray-400 hover:bg-gray-800 hover:text-gray-200"
            }`}
            onClick={() => {
              setActiveSequenceId(s.id);
              setActivePanel("editor");
            }}
          >
            <span className="text-gray-600 shrink-0">
              {s.molecule_type === "DNA" ? "🔵" : s.molecule_type === "RNA" ? "🟡" : "🟣"}
            </span>
            <div className="flex-1 min-w-0">
              <div className="truncate font-mono">{s.id}</div>
              <div className="text-gray-600">{s.length.toLocaleString()} bp</div>
            </div>
            <button
              onClick={async (e) => {
                e.stopPropagation();
                await seqApi.delete(s.id);
                removeSequence(s.id);
                qc.invalidateQueries({ queryKey: ["sequence", s.id] });
              }}
              className="opacity-0 group-hover:opacity-100 text-red-500 hover:text-red-400 text-sm"
            >
              ×
            </button>
          </div>
        ))}
      </div>

      {/* Status bar */}
      <div className="p-3 border-t border-gray-700 text-xs text-gray-500">
        <div className="flex justify-between">
          <span>{sequences.length} sequences</span>
          <span className="text-green-600">● API connected</span>
        </div>
      </div>
    </div>
  );
}
