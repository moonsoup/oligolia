import { useState } from "react";
import { useMutation } from "@tanstack/react-query";
import { api } from "../../api";
import type { AlignmentResult } from "../../types";

export function AlignmentViewer() {
  const [seq1, setSeq1] = useState("");
  const [seq2, setSeq2] = useState("");
  const [mode, setMode] = useState<"global" | "local">("global");
  const [msaSeqs, setMsaSeqs] = useState([{ id: "seq1", seq: "" }, { id: "seq2", seq: "" }]);

  const pairMut = useMutation({
    mutationFn: () => api.post<AlignmentResult>("/alignment/pairwise", { seq1, seq2, mode }).then((r) => r.data),
  });

  const msaMut = useMutation({
    mutationFn: () =>
      api.post<{ aligned: { id: string; aligned_seq: string }[]; consensus: string; identity_matrix: number[][] }>(
        "/alignment/multiple",
        { sequences: msaSeqs.filter((s) => s.seq), algorithm: "muscle" }
      ).then((r) => r.data),
  });

  return (
    <div className="p-4 space-y-4 h-full overflow-auto">
      {/* Pairwise */}
      <div className="bg-gray-800 rounded-lg p-4 space-y-3">
        <h3 className="font-semibold text-gray-200">Pairwise Alignment</h3>
        <div className="flex gap-2">
          <button
            onClick={() => setMode("global")}
            className={`px-3 py-1 rounded text-sm ${mode === "global" ? "bg-blue-700 text-white" : "bg-gray-700 text-gray-400"}`}
          >Global (Needleman-Wunsch)</button>
          <button
            onClick={() => setMode("local")}
            className={`px-3 py-1 rounded text-sm ${mode === "local" ? "bg-blue-700 text-white" : "bg-gray-700 text-gray-400"}`}
          >Local (Smith-Waterman)</button>
        </div>
        <textarea
          value={seq1}
          onChange={(e) => setSeq1(e.target.value)}
          placeholder="Sequence 1…"
          rows={2}
          className="w-full bg-gray-900 text-green-300 font-mono text-xs rounded p-2 resize-y"
        />
        <textarea
          value={seq2}
          onChange={(e) => setSeq2(e.target.value)}
          placeholder="Sequence 2…"
          rows={2}
          className="w-full bg-gray-900 text-green-300 font-mono text-xs rounded p-2 resize-y"
        />
        <button
          onClick={() => pairMut.mutate()}
          disabled={!seq1 || !seq2 || pairMut.isPending}
          className="bg-blue-600 hover:bg-blue-700 disabled:opacity-50 text-white px-4 py-2 rounded text-sm transition-colors"
        >
          {pairMut.isPending ? "Aligning…" : "Align"}
        </button>

        {pairMut.data && (
          <div className="space-y-2 text-xs">
            <div className="grid grid-cols-4 gap-3 text-center">
              {[
                ["Score", pairMut.data.score.toFixed(1)],
                ["Identity", pairMut.data.identity + "%"],
                ["Similarity", pairMut.data.similarity + "%"],
                ["Gaps", pairMut.data.gaps],
              ].map(([label, val]) => (
                <div key={label} className="bg-gray-900 rounded p-2">
                  <div className="text-gray-500">{label}</div>
                  <div className="text-white font-medium text-sm">{val}</div>
                </div>
              ))}
            </div>
            <div className="font-mono bg-gray-900 rounded p-2 overflow-x-auto">
              <div className="text-blue-300">{pairMut.data.aligned_seq1}</div>
              <div className="text-green-300">{pairMut.data.aligned_seq2}</div>
            </div>
          </div>
        )}
      </div>

      {/* MSA */}
      <div className="bg-gray-800 rounded-lg p-4 space-y-3">
        <div className="flex items-center justify-between">
          <h3 className="font-semibold text-gray-200">Multiple Sequence Alignment</h3>
          <button
            onClick={() => setMsaSeqs((s) => [...s, { id: `seq${s.length + 1}`, seq: "" }])}
            className="text-xs bg-gray-700 hover:bg-gray-600 text-white px-2 py-1 rounded"
          >
            + Add sequence
          </button>
        </div>
        {msaSeqs.map((s, i) => (
          <div key={i} className="flex gap-2">
            <input
              value={s.id}
              onChange={(e) => setMsaSeqs((seqs) => seqs.map((sq, j) => j === i ? { ...sq, id: e.target.value } : sq))}
              placeholder="ID"
              className="w-20 bg-gray-700 text-white rounded px-2 py-1 text-xs"
            />
            <input
              value={s.seq}
              onChange={(e) => setMsaSeqs((seqs) => seqs.map((sq, j) => j === i ? { ...sq, seq: e.target.value } : sq))}
              placeholder="Sequence (DNA/RNA/protein)…"
              className="flex-1 bg-gray-900 text-green-300 font-mono text-xs rounded px-2 py-1"
            />
          </div>
        ))}
        <button
          onClick={() => msaMut.mutate()}
          disabled={msaSeqs.filter((s) => s.seq).length < 2 || msaMut.isPending}
          className="bg-purple-600 hover:bg-purple-700 disabled:opacity-50 text-white px-4 py-2 rounded text-sm transition-colors"
        >
          {msaMut.isPending ? "Aligning…" : "Run MSA"}
        </button>

        {msaMut.data && (
          <div className="space-y-2 text-xs">
            <div className="font-mono bg-gray-900 rounded p-2 overflow-x-auto space-y-0.5">
              {msaMut.data.aligned.map((a) => (
                <div key={a.id} className="flex gap-2">
                  <span className="text-gray-400 w-16 shrink-0 truncate">{a.id}</span>
                  <span className="text-green-300 tracking-wider">{a.aligned_seq}</span>
                </div>
              ))}
              <div className="flex gap-2 border-t border-gray-700 pt-1">
                <span className="text-gray-500 w-16 shrink-0">consensus</span>
                <span className="text-yellow-300 tracking-wider">{msaMut.data.consensus}</span>
              </div>
            </div>
            <div>
              <p className="text-gray-500 mb-1">Pairwise identity matrix (%)</p>
              <div className="overflow-x-auto">
                <table className="text-center text-xs">
                  <thead>
                    <tr>
                      <th className="px-2 py-1" />
                      {msaMut.data.aligned.map((a) => (
                        <th key={a.id} className="px-2 py-1 text-gray-400">{a.id}</th>
                      ))}
                    </tr>
                  </thead>
                  <tbody>
                    {msaMut.data.identity_matrix.map((row, i) => (
                      <tr key={i}>
                        <td className="px-2 py-1 text-gray-400">{msaMut.data!.aligned[i].id}</td>
                        {row.map((val, j) => (
                          <td key={j} className={`px-2 py-1 ${val === 100 ? "text-green-400" : val > 80 ? "text-blue-400" : "text-gray-300"}`}>
                            {val}
                          </td>
                        ))}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
