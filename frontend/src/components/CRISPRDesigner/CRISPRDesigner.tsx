import { useState } from "react";
import { useMutation } from "@tanstack/react-query";
import { crispr } from "../../api";
import type { CasType, GuideRNA } from "../../types";
import { useAppStore } from "../../stores/appStore";
import { sequences } from "../../api";
import { useQuery } from "@tanstack/react-query";

const CAS_TYPES: { value: CasType; label: string; description: string }[] = [
  { value: "SpCas9", label: "SpCas9", description: "Most common — NGG PAM, 20 nt guide" },
  { value: "SpCas9-HF1", label: "SpCas9-HF1", description: "High-fidelity variant — fewer off-targets" },
  { value: "AsCas12a", label: "AsCas12a (Cpf1)", description: "TTTV PAM, 23 nt guide, staggered cuts" },
  { value: "LwaCas13a", label: "LwaCas13a", description: "RNA-targeting, no PAM required" },
];

function ScoreBar({ score }: { score: number }) {
  const pct = Math.round(score * 100);
  const color = score >= 0.6 ? "bg-green-500" : score >= 0.4 ? "bg-yellow-500" : "bg-red-500";
  return (
    <div className="flex items-center gap-2">
      <div className="flex-1 bg-gray-700 rounded-full h-1.5">
        <div className={`h-1.5 rounded-full ${color}`} style={{ width: `${pct}%` }} />
      </div>
      <span className="text-xs text-gray-400 w-8">{pct}%</span>
    </div>
  );
}

export function CRISPRDesigner() {
  const activeId = useAppStore((s) => s.activeSequenceId);
  const [casType, setCasType] = useState<CasType>("SpCas9");
  const [customTarget, setCustomTarget] = useState("");
  const [maxGuides, setMaxGuides] = useState(10);
  const [selected, setSelected] = useState<GuideRNA | null>(null);

  const { data: seq } = useQuery({
    queryKey: ["sequence", activeId],
    queryFn: () => sequences.get(activeId!),
    enabled: !!activeId,
  });

  const target = customTarget || (seq?.seq ?? "");

  const designMut = useMutation({
    mutationFn: () =>
      crispr.design({
        target_sequence: target,
        cas_type: casType,
        guide_length: 20,
        max_guides: maxGuides,
        check_off_targets: false,
      }),
  });

  return (
    <div className="p-4 space-y-4 h-full overflow-auto">
      <div className="bg-gray-800 rounded-lg p-4 space-y-3">
        <h3 className="font-semibold text-gray-200">Target Sequence</h3>
        {seq && (
          <div className="text-xs text-green-400 mb-1">
            Using loaded sequence: <span className="font-medium">{seq.name || seq.id}</span> ({seq.length} bp)
          </div>
        )}
        <textarea
          value={customTarget}
          onChange={(e) => setCustomTarget(e.target.value)}
          placeholder={seq ? "Leave blank to use loaded sequence, or paste a different target…" : "Paste target DNA sequence here…"}
          rows={4}
          className="w-full bg-gray-900 text-green-300 font-mono text-xs rounded p-2 resize-y"
        />
        <div className="text-xs text-gray-500">
          Target length: {target.length} bp
        </div>

        <div>
          <label className="text-sm text-gray-300 mb-1 block">Cas Enzyme</label>
          <div className="grid grid-cols-2 gap-2">
            {CAS_TYPES.map((ct) => (
              <button
                key={ct.value}
                onClick={() => setCasType(ct.value)}
                className={`text-left p-2 rounded border text-xs transition-colors ${
                  casType === ct.value
                    ? "border-green-500 bg-green-900/30 text-green-300"
                    : "border-gray-600 bg-gray-700 text-gray-400"
                }`}
              >
                <div className="font-medium">{ct.label}</div>
                <div className="text-gray-500 mt-0.5">{ct.description}</div>
              </button>
            ))}
          </div>
        </div>

        <div className="flex items-center gap-3">
          <label className="text-sm text-gray-300">Max guides:</label>
          <input
            type="number"
            value={maxGuides}
            onChange={(e) => setMaxGuides(parseInt(e.target.value))}
            min={1} max={50}
            className="w-20 bg-gray-700 text-white rounded px-2 py-1 text-sm"
          />
        </div>

        <button
          onClick={() => designMut.mutate()}
          disabled={!target || designMut.isPending}
          className="w-full bg-green-600 hover:bg-green-700 disabled:opacity-50 text-white rounded px-4 py-2 text-sm font-medium transition-colors"
        >
          {designMut.isPending ? "Designing guides…" : "Design Guide RNAs"}
        </button>
      </div>

      {designMut.isError && (
        <div className="text-red-400 text-sm bg-red-900/20 rounded p-2">
          {(designMut.error as Error).message}
        </div>
      )}

      {designMut.data && (
        <div className="space-y-2">
          <div className="flex justify-between items-center">
            <h3 className="font-semibold text-gray-200">
              Guide RNAs ({designMut.data.guides.length} of {designMut.data.total_candidates} candidates)
            </h3>
            <button
              onClick={() => {
                const tsv = designMut.data!.guides
                  .map((g) => `${g.sequence}\t${g.pam}\t${g.position}\t${g.strand}\t${g.gc_content}\t${g.on_target_score}`)
                  .join("\n");
                const blob = new Blob(["sequence\tpam\tposition\tstrand\tgc\tscore\n" + tsv], { type: "text/tab-separated-values" });
                const a = document.createElement("a");
                a.href = URL.createObjectURL(blob);
                a.download = "crispr_guides.tsv";
                a.click();
              }}
              className="text-xs bg-gray-700 hover:bg-gray-600 text-white px-3 py-1 rounded"
            >
              Export TSV
            </button>
          </div>

          {designMut.data.guides.map((g, i) => (
            <div
              key={i}
              onClick={() => setSelected(g === selected ? null : g)}
              className={`bg-gray-800 rounded-lg p-3 cursor-pointer transition-colors ${
                selected === g ? "border border-green-500" : "hover:bg-gray-750"
              }`}
            >
              <div className="flex items-center justify-between mb-1">
                <div className="font-mono text-sm text-green-300">
                  {g.sequence}
                  <span className="text-gray-500 ml-1 text-xs">|{g.pam}</span>
                </div>
                <span className={`text-xs px-1.5 py-0.5 rounded ${
                  g.strand === "+" ? "bg-blue-900 text-blue-300" : "bg-orange-900 text-orange-300"
                }`}>
                  {g.strand} strand
                </span>
              </div>
              <div className="grid grid-cols-3 gap-2 text-xs text-gray-400 mb-1">
                <span>Pos: {g.position}</span>
                <span>GC: {g.gc_content}%</span>
              </div>
              {g.on_target_score != null && (
                <div>
                  <div className="text-xs text-gray-500 mb-0.5">On-target score</div>
                  <ScoreBar score={g.on_target_score} />
                </div>
              )}
              {selected === g && (
                <div className="mt-2 pt-2 border-t border-gray-700 space-y-1">
                  <div className="text-xs text-gray-400">
                    <strong>Full guide + PAM context:</strong>
                    <div className="font-mono text-green-300 mt-1">5'–{g.sequence}{g.pam}–3'</div>
                  </div>
                  <button
                    onClick={(e) => { e.stopPropagation(); navigator.clipboard.writeText(g.sequence); }}
                    className="text-xs bg-gray-700 hover:bg-gray-600 text-white px-2 py-1 rounded"
                  >
                    Copy guide
                  </button>
                </div>
              )}
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
