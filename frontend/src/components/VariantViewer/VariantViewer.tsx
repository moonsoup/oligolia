import { useState } from "react";
import { useAppStore } from "../../stores/appStore";
import type { Variant } from "../../types";

const SIG_COLORS: Record<string, string> = {
  Pathogenic: "text-red-400",
  "Likely pathogenic": "text-orange-400",
  "Uncertain significance": "text-yellow-400",
  "Likely benign": "text-blue-400",
  Benign: "text-green-400",
};

function VariantRow({ v, onClick }: { v: Variant; onClick: () => void }) {
  return (
    <tr
      className="border-b border-gray-700 hover:bg-gray-750 cursor-pointer text-xs"
      onClick={onClick}
    >
      <td className="px-2 py-1.5 text-gray-300">{v.chrom}:{v.pos}</td>
      <td className="px-2 py-1.5 font-mono text-red-300">{v.ref}</td>
      <td className="px-2 py-1.5 font-mono text-green-300">{v.alt.join(",")}</td>
      <td className="px-2 py-1.5 text-blue-300">{v.gene || "—"}</td>
      <td className="px-2 py-1.5">
        <span className={`px-1.5 py-0.5 rounded text-xs ${
          v.variant_type === "SNP" ? "bg-blue-900 text-blue-300" :
          v.variant_type === "DEL" ? "bg-red-900 text-red-300" :
          v.variant_type === "INS" ? "bg-green-900 text-green-300" :
          "bg-gray-700 text-gray-300"
        }`}>{v.variant_type}</span>
      </td>
      <td className={`px-2 py-1.5 ${SIG_COLORS[v.clinical_significance ?? ""] ?? "text-gray-400"}`}>
        {v.clinical_significance ?? "—"}
      </td>
      <td className="px-2 py-1.5 text-gray-400">
        {v.gnomad_af != null ? v.gnomad_af.toExponential(2) : "—"}
      </td>
    </tr>
  );
}

export function VariantViewer() {
  const { variants } = useAppStore((s) => ({ variants: s.variants }));
  const [selected, setSelected] = useState<Variant | null>(null);
  const [filter, setFilter] = useState("");

  const filtered = variants.filter(
    (v) =>
      !filter ||
      v.chrom.includes(filter) ||
      v.gene.toLowerCase().includes(filter.toLowerCase()) ||
      v.ref.includes(filter.toUpperCase()) ||
      v.alt.some((a) => a.includes(filter.toUpperCase()))
  );

  if (!variants.length) {
    return (
      <div className="flex flex-col items-center justify-center h-full text-gray-400 space-y-2">
        <p>No variants loaded.</p>
        <p className="text-xs">Upload a VCF file via the File Manager to get started.</p>
      </div>
    );
  }

  return (
    <div className="flex flex-col h-full">
      <div className="p-4 border-b border-gray-700 flex gap-2">
        <input
          value={filter}
          onChange={(e) => setFilter(e.target.value)}
          placeholder="Filter by chromosome, gene, or allele…"
          className="flex-1 bg-gray-700 text-white rounded px-3 py-2 text-sm"
        />
        <span className="text-gray-400 text-sm self-center">
          {filtered.length}/{variants.length}
        </span>
      </div>

      <div className="flex-1 overflow-auto">
        <table className="w-full border-collapse">
          <thead className="sticky top-0 bg-gray-900">
            <tr className="text-left text-xs text-gray-400 border-b border-gray-700">
              <th className="px-2 py-2">Position</th>
              <th className="px-2 py-2">Ref</th>
              <th className="px-2 py-2">Alt</th>
              <th className="px-2 py-2">Gene</th>
              <th className="px-2 py-2">Type</th>
              <th className="px-2 py-2">Clinical Sig.</th>
              <th className="px-2 py-2">gnomAD AF</th>
            </tr>
          </thead>
          <tbody>
            {filtered.map((v, i) => (
              <VariantRow
                key={i}
                v={v}
                onClick={() => setSelected(v === selected ? null : v)}
              />
            ))}
          </tbody>
        </table>
      </div>

      {selected && (
        <div className="border-t border-gray-700 bg-gray-800 p-4 max-h-48 overflow-auto">
          <div className="flex justify-between items-start">
            <h3 className="text-sm font-semibold text-gray-200">Variant Detail</h3>
            <button
              onClick={() => setSelected(null)}
              className="text-gray-500 hover:text-gray-300 text-lg"
            >
              ×
            </button>
          </div>
          <div className="grid grid-cols-2 gap-x-8 gap-y-1 mt-2 text-xs text-gray-300">
            <div><span className="text-gray-500">ID:</span> {selected.id}</div>
            <div><span className="text-gray-500">Gene:</span> {selected.gene || "—"}</div>
            <div><span className="text-gray-500">Position:</span> {selected.chrom}:{selected.pos}</div>
            <div><span className="text-gray-500">Ref→Alt:</span> {selected.ref}→{selected.alt.join(",")}</div>
            <div><span className="text-gray-500">QUAL:</span> {selected.qual ?? "—"}</div>
            <div><span className="text-gray-500">Filter:</span> {selected.filter.join(";") || "PASS"}</div>
            {Object.entries(selected.info).map(([k, v]) => (
              <div key={k}><span className="text-gray-500">{k}:</span> {String(v)}</div>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
