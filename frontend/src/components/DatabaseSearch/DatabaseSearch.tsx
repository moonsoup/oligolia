import { useState } from "react";
import { useMutation } from "@tanstack/react-query";
import { databases, sequences as seqApi } from "../../api";
import type { Database, SearchResult } from "../../types";
import { useAppStore } from "../../stores/appStore";

const ALL_DATABASES: { value: Database; label: string }[] = [
  { value: "ncbi_gene", label: "NCBI Gene" },
  { value: "ncbi_nucleotide", label: "NCBI Nucleotide" },
  { value: "ensembl", label: "Ensembl" },
  { value: "uniprot", label: "UniProt" },
  { value: "kegg", label: "KEGG" },
];

const DB_COLORS: Record<Database, string> = {
  ncbi_gene: "bg-blue-800 text-blue-200",
  ncbi_nucleotide: "bg-blue-700 text-blue-100",
  ncbi_protein: "bg-indigo-800 text-indigo-200",
  ncbi_snp: "bg-purple-800 text-purple-200",
  ncbi_clinvar: "bg-red-800 text-red-200",
  ensembl: "bg-green-800 text-green-200",
  uniprot: "bg-yellow-800 text-yellow-200",
  kegg: "bg-orange-800 text-orange-200",
  reactome: "bg-pink-800 text-pink-200",
  gnomad: "bg-teal-800 text-teal-200",
  string: "bg-cyan-800 text-cyan-200",
  pdb: "bg-gray-700 text-gray-200",
};

export function DatabaseSearch() {
  const [query, setQuery] = useState("");
  const [species, setSpecies] = useState("homo sapiens");
  const [selectedDBs, setSelectedDBs] = useState<Database[]>(["ncbi_gene", "ensembl"]);
  const addSeq = useAppStore((s) => s.addSequence);
  const setActiveId = useAppStore((s) => s.setActiveSequenceId);
  const setPanel = useAppStore((s) => s.setActivePanel);

  const searchMut = useMutation({
    mutationFn: () =>
      databases.search({ query, databases: selectedDBs, species, max_results: 20 }),
  });

  const fetchMut = useMutation({
    mutationFn: async (result: SearchResult) => {
      if (result.database.startsWith("ncbi")) {
        const data = await databases.fetchNCBISequence(result.accession, "fasta");
        return { id: result.accession, name: result.name, seq: data.content, description: result.description };
      }
      if (result.database === "ensembl") {
        const data = await databases.fetchEnsemblSequence(result.accession, "cdna");
        return { id: result.accession, name: result.name, seq: data.sequence, description: data.desc };
      }
      throw new Error(`Sequence fetch not available for ${result.database}`);
    },
    onSuccess: async (data) => {
      const seq = {
        id: data.id,
        name: data.name,
        description: data.description,
        seq: data.seq.replace(/^>.*\n/, "").replace(/\n/g, ""),
        molecule_type: "DNA" as const,
        annotations: [],
        source_db: "",
        accession: data.id,
        length: 0,
      };
      await seqApi.add(seq);
      addSeq({ ...seq, length: seq.seq.length });
      setActiveId(seq.id);
      setPanel("editor");
    },
  });

  const toggleDB = (db: Database) => {
    setSelectedDBs((prev) =>
      prev.includes(db) ? prev.filter((d) => d !== db) : [...prev, db]
    );
  };

  return (
    <div className="p-4 space-y-4 h-full overflow-auto">
      <div className="space-y-2">
        <div className="flex gap-2">
          <input
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && searchMut.mutate()}
            placeholder="Search gene (e.g. BRCA2, TP53, CFTR)…"
            className="flex-1 bg-gray-700 text-white rounded px-3 py-2 text-sm"
          />
          <button
            onClick={() => searchMut.mutate()}
            disabled={!query || searchMut.isPending}
            className="bg-blue-600 hover:bg-blue-700 disabled:opacity-50 text-white px-4 py-2 rounded text-sm transition-colors"
          >
            {searchMut.isPending ? "Searching…" : "Search"}
          </button>
        </div>
        <input
          value={species}
          onChange={(e) => setSpecies(e.target.value)}
          placeholder="Species (e.g. homo sapiens)"
          className="w-full bg-gray-700 text-white rounded px-3 py-2 text-sm"
        />
        <div className="flex flex-wrap gap-2">
          {ALL_DATABASES.map((db) => (
            <button
              key={db.value}
              onClick={() => toggleDB(db.value)}
              className={`text-xs px-2 py-1 rounded border transition-colors ${
                selectedDBs.includes(db.value)
                  ? "border-blue-500 bg-blue-900 text-blue-200"
                  : "border-gray-600 bg-gray-800 text-gray-400"
              }`}
            >
              {db.label}
            </button>
          ))}
        </div>
      </div>

      {searchMut.isError && (
        <div className="text-red-400 text-sm bg-red-900/20 rounded p-2">
          Search failed. Check connection and try again.
        </div>
      )}

      {searchMut.data && (
        <div className="space-y-2">
          <p className="text-gray-400 text-xs">
            {searchMut.data.total} results from {searchMut.data.databases_searched.length} databases
          </p>
          {Object.entries(searchMut.data.errors).length > 0 && (
            <div className="text-yellow-400 text-xs">
              Errors: {Object.entries(searchMut.data.errors).map(([db, err]) => `${db}: ${err}`).join("; ")}
            </div>
          )}
          {searchMut.data.results.map((r) => (
            <div key={`${r.database}-${r.id}`} className="bg-gray-800 rounded-lg p-3 space-y-1">
              <div className="flex items-start justify-between gap-2">
                <div>
                  <span className="text-green-400 font-medium text-sm">{r.name}</span>
                  <span className={`ml-2 text-xs px-1.5 py-0.5 rounded ${DB_COLORS[r.database]}`}>
                    {r.database}
                  </span>
                </div>
                <button
                  onClick={() => fetchMut.mutate(r)}
                  disabled={fetchMut.isPending}
                  className="text-xs bg-green-800 hover:bg-green-700 text-green-200 px-2 py-1 rounded transition-colors shrink-0"
                >
                  Load
                </button>
              </div>
              <p className="text-gray-400 text-xs line-clamp-2">{r.description}</p>
              {r.organism && <p className="text-gray-500 text-xs italic">{r.organism}</p>}
              {r.url && (
                <a href={r.url} target="_blank" rel="noreferrer"
                  className="text-blue-400 text-xs hover:underline">{r.url}</a>
              )}
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
