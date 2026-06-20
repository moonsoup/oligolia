import { useState } from "react";
import { useMutation } from "@tanstack/react-query";
import { api } from "../../api";
import type { PrimerPair, RestrictionSite } from "../../types";
import { useAppStore } from "../../stores/appStore";
import { sequences } from "../../api";
import { useQuery } from "@tanstack/react-query";

export function PrimerDesigner() {
  const activeId = useAppStore((s) => s.activeSequenceId);
  const { data: seq } = useQuery({
    queryKey: ["sequence", activeId],
    queryFn: () => sequences.get(activeId!),
    enabled: !!activeId,
  });

  const [customTemplate, setCustomTemplate] = useState("");
  const [productMin, setProductMin] = useState(100);
  const [productMax, setProductMax] = useState(600);
  const [tmMin, setTmMin] = useState(55);
  const [tmMax, setTmMax] = useState(65);

  const template = customTemplate || (seq?.seq ?? "");

  const primerMut = useMutation({
    mutationFn: () =>
      api.post<PrimerPair[]>("/primers/design", {
        template,
        product_min: productMin,
        product_max: productMax,
        tm_min: tmMin,
        tm_max: tmMax,
        max_pairs: 5,
      }).then((r) => r.data),
  });

  const restrictionMut = useMutation({
    mutationFn: () =>
      api.post<RestrictionSite[]>("/primers/restriction_sites", { template }).then((r) => r.data),
  });

  return (
    <div className="p-4 space-y-4 h-full overflow-auto">
      <div className="bg-gray-800 rounded-lg p-4 space-y-3">
        <h3 className="font-semibold text-gray-200">Template</h3>
        {seq && (
          <div className="text-xs text-green-400">
            Using: <span className="font-medium">{seq.name || seq.id}</span>
          </div>
        )}
        <textarea
          value={customTemplate}
          onChange={(e) => setCustomTemplate(e.target.value)}
          placeholder={seq ? "Leave blank to use loaded sequence…" : "Paste template sequence…"}
          rows={3}
          className="w-full bg-gray-900 text-green-300 font-mono text-xs rounded p-2 resize-y"
        />
      </div>

      {/* PCR Primer Design */}
      <div className="bg-gray-800 rounded-lg p-4 space-y-3">
        <h3 className="font-semibold text-gray-200">PCR Primer Design</h3>
        <div className="grid grid-cols-2 gap-2 text-xs">
          <div>
            <label className="text-gray-400">Product min (bp)</label>
            <input type="number" value={productMin} onChange={(e) => setProductMin(+e.target.value)}
              className="w-full bg-gray-700 text-white rounded px-2 py-1 mt-1" />
          </div>
          <div>
            <label className="text-gray-400">Product max (bp)</label>
            <input type="number" value={productMax} onChange={(e) => setProductMax(+e.target.value)}
              className="w-full bg-gray-700 text-white rounded px-2 py-1 mt-1" />
          </div>
          <div>
            <label className="text-gray-400">Tm min (°C)</label>
            <input type="number" value={tmMin} onChange={(e) => setTmMin(+e.target.value)}
              className="w-full bg-gray-700 text-white rounded px-2 py-1 mt-1" />
          </div>
          <div>
            <label className="text-gray-400">Tm max (°C)</label>
            <input type="number" value={tmMax} onChange={(e) => setTmMax(+e.target.value)}
              className="w-full bg-gray-700 text-white rounded px-2 py-1 mt-1" />
          </div>
        </div>
        <button
          onClick={() => primerMut.mutate()}
          disabled={!template || primerMut.isPending}
          className="w-full bg-orange-600 hover:bg-orange-700 disabled:opacity-50 text-white py-2 rounded text-sm transition-colors"
        >
          {primerMut.isPending ? "Designing…" : "Design Primers"}
        </button>

        {primerMut.data?.map((pair, i) => (
          <div key={i} className="bg-gray-900 rounded p-3 text-xs space-y-1.5">
            <div className="flex justify-between text-gray-400">
              <span>Pair {i + 1}</span>
              <span>Product: <strong className="text-white">{pair.product_size} bp</strong></span>
            </div>
            <div className="font-mono">
              <span className="text-blue-300">F: {pair.forward.sequence}</span>
              <span className="text-gray-500 ml-2">Tm {pair.forward.tm}°C · GC {pair.forward.gc_content}%</span>
            </div>
            <div className="font-mono">
              <span className="text-orange-300">R: {pair.reverse.sequence}</span>
              <span className="text-gray-500 ml-2">Tm {pair.reverse.tm}°C · GC {pair.reverse.gc_content}%</span>
            </div>
          </div>
        ))}
      </div>

      {/* Restriction enzyme sites */}
      <div className="bg-gray-800 rounded-lg p-4 space-y-3">
        <h3 className="font-semibold text-gray-200">Restriction Enzyme Analysis</h3>
        <button
          onClick={() => restrictionMut.mutate()}
          disabled={!template || restrictionMut.isPending}
          className="w-full bg-teal-600 hover:bg-teal-700 disabled:opacity-50 text-white py-2 rounded text-sm transition-colors"
        >
          {restrictionMut.isPending ? "Analysing…" : "Find Cut Sites"}
        </button>
        {restrictionMut.data && (
          <div className="space-y-1">
            {restrictionMut.data.length === 0 ? (
              <p className="text-gray-500 text-xs">No restriction sites found.</p>
            ) : (
              restrictionMut.data.map((s) => (
                <div key={s.enzyme} className="flex items-center justify-between text-xs">
                  <span className="text-green-400 font-medium w-20">{s.enzyme}</span>
                  <span className="font-mono text-gray-400">{s.cut_pattern}</span>
                  <span className="text-blue-300">{s.count}×</span>
                  <span className="text-gray-500 truncate max-w-24">
                    {s.positions.slice(0, 4).join(", ")}{s.positions.length > 4 ? "…" : ""}
                  </span>
                </div>
              ))
            )}
          </div>
        )}
      </div>
    </div>
  );
}
