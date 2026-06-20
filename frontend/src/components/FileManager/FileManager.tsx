import { useRef, useState } from "react";
import { api, sequences as seqApi } from "../../api";
import type { Sequence, Variant } from "../../types";
import { useAppStore } from "../../stores/appStore";

type Tab = "import" | "export" | "share";

export function FileManager() {
  const [tab, setTab] = useState<Tab>("import");
  const [importing, setImporting] = useState(false);
  const [importError, setImportError] = useState<string | null>(null);
  const [importedSeqs, setImportedSeqs] = useState<Sequence[]>([]);
  const fileRef = useRef<HTMLInputElement>(null);
  const { sequences, variants, addSequence } = useAppStore((s) => ({
    sequences: s.sequences,
    variants: s.variants,
    addSequence: s.addSequence,
  }));
  const [emailTarget, setEmailTarget] = useState("");
  const [emailStatus, setEmailStatus] = useState("");

  const handleFileImport = async (file: File) => {
    setImporting(true);
    setImportError(null);
    try {
      const ext = file.name.split(".").pop()?.toLowerCase();
      const formData = new FormData();
      formData.append("file", file);

      if (ext === "vcf") {
        const r = await api.post<Variant[]>("/files/parse/vcf", formData, {
          headers: { "Content-Type": "multipart/form-data" },
        });
        useAppStore.getState().setVariants(r.data);
        setImportedSeqs([]);
      } else if (ext === "gff3" || ext === "gff" || ext === "gtf") {
        await api.post("/files/parse/gff", formData, {
          headers: { "Content-Type": "multipart/form-data" },
        });
      } else {
        const r = await api.post<Sequence[]>("/files/parse/sequence", formData, {
          headers: { "Content-Type": "multipart/form-data" },
        });
        for (const seq of r.data) {
          await seqApi.add(seq);
          addSequence(seq);
        }
        setImportedSeqs(r.data);
      }
    } catch (e: unknown) {
      setImportError(e instanceof Error ? e.message : "Import failed");
    } finally {
      setImporting(false);
    }
  };

  const downloadFasta = async () => {
    if (!sequences.length) return;
    const r = await api.post("/files/download/fasta", sequences, { responseType: "blob" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(r.data as Blob);
    a.download = "sequences.fa";
    a.click();
  };

  const downloadGenbank = async () => {
    if (!sequences.length) return;
    const r = await api.post("/files/download/genbank", sequences, { responseType: "blob" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(r.data as Blob);
    a.download = "sequences.gb";
    a.click();
  };

  const downloadVcf = async () => {
    if (!variants.length) return;
    const r = await api.post("/files/download/vcf", variants, { responseType: "blob" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(r.data as Blob);
    a.download = "variants.vcf";
    a.click();
  };

  return (
    <div className="p-4 space-y-4 h-full overflow-auto">
      <div className="flex gap-1 bg-gray-900 rounded-lg p-1">
        {(["import", "export", "share"] as Tab[]).map((t) => (
          <button
            key={t}
            onClick={() => setTab(t)}
            className={`flex-1 py-1.5 rounded text-sm font-medium transition-colors ${
              tab === t ? "bg-gray-700 text-white" : "text-gray-400 hover:text-gray-200"
            }`}
          >
            {t.charAt(0).toUpperCase() + t.slice(1)}
          </button>
        ))}
      </div>

      {tab === "import" && (
        <div className="space-y-4">
          <div
            className="border-2 border-dashed border-gray-600 rounded-lg p-8 text-center cursor-pointer hover:border-green-500 transition-colors"
            onClick={() => fileRef.current?.click()}
            onDragOver={(e) => e.preventDefault()}
            onDrop={(e) => {
              e.preventDefault();
              const file = e.dataTransfer.files[0];
              if (file) handleFileImport(file);
            }}
          >
            <input
              ref={fileRef}
              type="file"
              className="hidden"
              accept=".fasta,.fa,.fna,.faa,.fastq,.fq,.gb,.gbk,.embl,.vcf,.gff,.gff3,.gtf"
              onChange={(e) => e.target.files?.[0] && handleFileImport(e.target.files[0])}
            />
            {importing ? (
              <p className="text-blue-400">Importing…</p>
            ) : (
              <>
                <p className="text-gray-300">Drop file here or click to browse</p>
                <p className="text-gray-500 text-xs mt-1">
                  FASTA, FASTQ, GenBank, EMBL, VCF, GFF3/GTF
                </p>
              </>
            )}
          </div>

          {importError && <p className="text-red-400 text-sm">{importError}</p>}

          {importedSeqs.length > 0 && (
            <div className="bg-gray-800 rounded-lg p-3">
              <p className="text-green-400 text-sm font-medium mb-2">
                Imported {importedSeqs.length} sequence{importedSeqs.length > 1 ? "s" : ""}
              </p>
              {importedSeqs.map((s) => (
                <div key={s.id} className="text-xs text-gray-400 flex justify-between">
                  <span className="font-mono">{s.id}</span>
                  <span>{s.length.toLocaleString()} bp</span>
                </div>
              ))}
            </div>
          )}
        </div>
      )}

      {tab === "export" && (
        <div className="space-y-3">
          <p className="text-gray-400 text-sm">
            {sequences.length} sequence{sequences.length !== 1 ? "s" : ""} loaded •{" "}
            {variants.length} variant{variants.length !== 1 ? "s" : ""} loaded
          </p>

          <div className="space-y-2">
            <button
              onClick={downloadFasta}
              disabled={!sequences.length}
              className="w-full bg-gray-700 hover:bg-gray-600 disabled:opacity-40 text-white px-4 py-3 rounded-lg text-sm text-left transition-colors"
            >
              <div className="font-medium">Download FASTA (.fa)</div>
              <div className="text-gray-400 text-xs mt-0.5">All loaded sequences</div>
            </button>

            <button
              onClick={downloadGenbank}
              disabled={!sequences.length}
              className="w-full bg-gray-700 hover:bg-gray-600 disabled:opacity-40 text-white px-4 py-3 rounded-lg text-sm text-left transition-colors"
            >
              <div className="font-medium">Download GenBank (.gb)</div>
              <div className="text-gray-400 text-xs mt-0.5">All loaded sequences with annotations</div>
            </button>

            <button
              onClick={downloadVcf}
              disabled={!variants.length}
              className="w-full bg-gray-700 hover:bg-gray-600 disabled:opacity-40 text-white px-4 py-3 rounded-lg text-sm text-left transition-colors"
            >
              <div className="font-medium">Download VCF (.vcf)</div>
              <div className="text-gray-400 text-xs mt-0.5">All loaded variants</div>
            </button>
          </div>
        </div>
      )}

      {tab === "share" && (
        <div className="space-y-4">
          <div className="bg-gray-800 rounded-lg p-4 space-y-3">
            <h3 className="font-semibold text-gray-200">Email Sequences</h3>
            <input
              value={emailTarget}
              onChange={(e) => setEmailTarget(e.target.value)}
              placeholder="Recipient email address"
              type="email"
              className="w-full bg-gray-700 text-white rounded px-3 py-2 text-sm"
            />
            <button
              disabled={!emailTarget || !sequences.length}
              onClick={async () => {
                try {
                  await api.post("/share/email", {
                    to: emailTarget,
                    sequences,
                  });
                  setEmailStatus("Email sent successfully!");
                } catch {
                  setEmailStatus("Email failed — configure SMTP in settings.");
                }
              }}
              className="w-full bg-blue-600 hover:bg-blue-700 disabled:opacity-40 text-white py-2 rounded text-sm transition-colors"
            >
              Send FASTA by Email
            </button>
            {emailStatus && <p className="text-sm text-green-400">{emailStatus}</p>}
          </div>

          <div className="bg-gray-800 rounded-lg p-4 space-y-2">
            <h3 className="font-semibold text-gray-200">Copy to Clipboard</h3>
            <button
              disabled={!sequences.length}
              onClick={() => {
                const fasta = sequences.map((s) => `>${s.id}\n${s.seq}`).join("\n");
                navigator.clipboard.writeText(fasta);
              }}
              className="w-full bg-gray-700 hover:bg-gray-600 disabled:opacity-40 text-white py-2 rounded text-sm transition-colors"
            >
              Copy all sequences as FASTA
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
