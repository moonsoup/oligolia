import { useState } from "react";
import { useMutation, useQuery } from "@tanstack/react-query";
import { sequences } from "../../api";
import type { Sequence, SequenceEditRequest } from "../../types";
import { useAppStore } from "../../stores/appStore";

const OPERATIONS = [
  { value: "reverse_complement", label: "Reverse Complement" },
  { value: "complement", label: "Complement" },
  { value: "translate", label: "Translate (DNA→Protein)" },
  { value: "transcribe", label: "Transcribe (DNA→RNA)" },
  { value: "back_transcribe", label: "Back-Transcribe (RNA→DNA)" },
  { value: "insert", label: "Insert Bases" },
  { value: "delete", label: "Delete Region" },
  { value: "replace", label: "Replace Region" },
];

export function SequenceEditor() {
  const activeId = useAppStore((s) => s.activeSequenceId);
  const [op, setOp] = useState<string>("reverse_complement");
  const [position, setPosition] = useState("");
  const [endPosition, setEndPosition] = useState("");
  const [insertSeq, setInsertSeq] = useState("");
  const [replacement, setReplacement] = useState("");
  const [result, setResult] = useState<string | null>(null);
  const [motif, setMotif] = useState("");
  const [motifResult, setMotifResult] = useState<{ count: number; positions: number[] } | null>(null);

  const { data: seq } = useQuery({
    queryKey: ["sequence", activeId],
    queryFn: () => sequences.get(activeId!),
    enabled: !!activeId,
  });

  const { data: gcData } = useQuery({
    queryKey: ["gc", activeId],
    queryFn: () => sequences.gcContent(activeId!),
    enabled: !!activeId,
  });

  const editMut = useMutation({
    mutationFn: (req: SequenceEditRequest) => sequences.edit(activeId!, req),
    onSuccess: (data) => setResult(data.result_seq),
  });

  const motifMut = useMutation({
    mutationFn: (m: string) => sequences.findMotif(activeId!, m),
    onSuccess: (data) => setMotifResult(data),
  });

  if (!activeId || !seq) {
    return (
      <div className="flex items-center justify-center h-full text-gray-400">
        <p>Select a sequence from the panel to edit it.</p>
      </div>
    );
  }

  const buildRequest = (): SequenceEditRequest => {
    const base: SequenceEditRequest = { operation: op as SequenceEditRequest["operation"] };
    if (op === "insert") {
      base.position = parseInt(position);
      base.insert_seq = insertSeq;
    }
    if (op === "delete") {
      base.position = parseInt(position);
      base.end_position = parseInt(endPosition);
    }
    if (op === "replace") {
      base.position = parseInt(position);
      base.end_position = parseInt(endPosition);
      base.replacement = replacement;
    }
    return base;
  };

  return (
    <div className="p-4 h-full overflow-auto space-y-4">
      {/* Sequence info */}
      <div className="bg-gray-800 rounded-lg p-4">
        <div className="flex justify-between items-start mb-2">
          <div>
            <h2 className="text-lg font-bold text-green-400">{seq.name || seq.id}</h2>
            <p className="text-gray-400 text-sm">{seq.description}</p>
          </div>
          <div className="text-right text-sm text-gray-400">
            <div>{seq.molecule_type}</div>
            <div>{seq.length.toLocaleString()} bp</div>
            {gcData && <div>GC: {gcData.gc_content}%</div>}
          </div>
        </div>
        <div className="font-mono text-xs text-green-300 bg-gray-900 rounded p-2 overflow-x-auto whitespace-nowrap">
          {seq.seq.length > 300 ? seq.seq.slice(0, 300) + "…" : seq.seq}
        </div>
      </div>

      {/* Annotations */}
      {seq.annotations.length > 0 && (
        <div className="bg-gray-800 rounded-lg p-4">
          <h3 className="font-semibold text-gray-200 mb-2">Annotations ({seq.annotations.length})</h3>
          <div className="space-y-1 max-h-32 overflow-y-auto">
            {seq.annotations.map((ann, i) => (
              <div key={i} className="text-xs flex gap-2 text-gray-300">
                <span className="text-blue-400 font-medium">{ann.feature_type}</span>
                <span>{ann.start}–{ann.end}</span>
                <span className="text-gray-500">{ann.strand}</span>
                {typeof ann.qualifiers.gene === "string" && (
                  <span className="text-yellow-400">{ann.qualifiers.gene}</span>
                )}
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Edit operations */}
      <div className="bg-gray-800 rounded-lg p-4 space-y-3">
        <h3 className="font-semibold text-gray-200">Edit Operation</h3>
        <select
          value={op}
          onChange={(e) => setOp(e.target.value)}
          className="w-full bg-gray-700 text-white rounded px-3 py-2 text-sm"
        >
          {OPERATIONS.map((o) => (
            <option key={o.value} value={o.value}>{o.label}</option>
          ))}
        </select>

        {["insert", "delete", "replace"].includes(op) && (
          <div className="grid grid-cols-2 gap-2">
            <input
              placeholder="Start position"
              value={position}
              onChange={(e) => setPosition(e.target.value)}
              className="bg-gray-700 text-white rounded px-3 py-2 text-sm"
            />
            {["delete", "replace"].includes(op) && (
              <input
                placeholder="End position"
                value={endPosition}
                onChange={(e) => setEndPosition(e.target.value)}
                className="bg-gray-700 text-white rounded px-3 py-2 text-sm"
              />
            )}
          </div>
        )}

        {op === "insert" && (
          <input
            placeholder="Sequence to insert (e.g. ATCG)"
            value={insertSeq}
            onChange={(e) => setInsertSeq(e.target.value)}
            className="w-full bg-gray-700 text-white rounded px-3 py-2 text-sm font-mono"
          />
        )}

        {op === "replace" && (
          <input
            placeholder="Replacement sequence"
            value={replacement}
            onChange={(e) => setReplacement(e.target.value)}
            className="w-full bg-gray-700 text-white rounded px-3 py-2 text-sm font-mono"
          />
        )}

        <button
          onClick={() => editMut.mutate(buildRequest())}
          disabled={editMut.isPending}
          className="w-full bg-green-600 hover:bg-green-700 disabled:opacity-50 text-white rounded px-4 py-2 text-sm font-medium transition-colors"
        >
          {editMut.isPending ? "Processing…" : "Apply Operation"}
        </button>

        {editMut.isError && (
          <p className="text-red-400 text-sm">{String((editMut.error as Error).message)}</p>
        )}
      </div>

      {/* Result */}
      {result && (
        <div className="bg-gray-800 rounded-lg p-4">
          <h3 className="font-semibold text-gray-200 mb-2">Result</h3>
          <div className="font-mono text-xs text-green-300 bg-gray-900 rounded p-2 overflow-x-auto whitespace-pre-wrap break-all">
            {result}
          </div>
          <div className="flex gap-2 mt-2">
            <button
              onClick={() => navigator.clipboard.writeText(result)}
              className="text-xs bg-gray-700 hover:bg-gray-600 text-white px-3 py-1 rounded transition-colors"
            >
              Copy
            </button>
            <button
              onClick={() => {
                const blob = new Blob([`>${activeId}_edited\n${result}`], { type: "text/plain" });
                const a = document.createElement("a");
                a.href = URL.createObjectURL(blob);
                a.download = `${activeId}_edited.fasta`;
                a.click();
              }}
              className="text-xs bg-blue-700 hover:bg-blue-600 text-white px-3 py-1 rounded transition-colors"
            >
              Download FASTA
            </button>
          </div>
        </div>
      )}

      {/* Motif search */}
      <div className="bg-gray-800 rounded-lg p-4 space-y-2">
        <h3 className="font-semibold text-gray-200">Find Motif (IUPAC)</h3>
        <div className="flex gap-2">
          <input
            placeholder="e.g. GAATTC or RRYYY"
            value={motif}
            onChange={(e) => setMotif(e.target.value)}
            className="flex-1 bg-gray-700 text-white rounded px-3 py-2 text-sm font-mono"
          />
          <button
            onClick={() => motifMut.mutate(motif)}
            disabled={!motif || motifMut.isPending}
            className="bg-purple-600 hover:bg-purple-700 disabled:opacity-50 text-white px-4 py-2 rounded text-sm transition-colors"
          >
            Search
          </button>
        </div>
        {motifResult && (
          <div className="text-sm text-gray-300">
            <span className="text-green-400 font-medium">{motifResult.count}</span> occurrences found
            {motifResult.positions.length > 0 && (
              <span className="text-gray-500 ml-2">
                at positions: {motifResult.positions.slice(0, 10).join(", ")}
                {motifResult.positions.length > 10 && "…"}
              </span>
            )}
          </div>
        )}
      </div>
    </div>
  );
}
