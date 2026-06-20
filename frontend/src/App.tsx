import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppStore } from "./stores/appStore";
import { Sidebar } from "./components/Sidebar/Sidebar";
import { SequenceEditor } from "./components/SequenceEditor/SequenceEditor";
import { DatabaseSearch } from "./components/DatabaseSearch/DatabaseSearch";
import { CRISPRDesigner } from "./components/CRISPRDesigner/CRISPRDesigner";
import { VariantViewer } from "./components/VariantViewer/VariantViewer";
import { AlignmentViewer } from "./components/AlignmentViewer/AlignmentViewer";
import { FileManager } from "./components/FileManager/FileManager";
import { PrimerDesigner } from "./components/PrimerDesigner/PrimerDesigner";

const queryClient = new QueryClient({
  defaultOptions: { queries: { retry: 1, staleTime: 30000 } },
});

function MainPanel() {
  const activePanel = useAppStore((s) => s.activePanel);

  switch (activePanel) {
    case "editor": return <SequenceEditor />;
    case "search": return <DatabaseSearch />;
    case "crispr": return <CRISPRDesigner />;
    case "variants": return <VariantViewer />;
    case "alignment": return <AlignmentViewer />;
    case "primers": return <PrimerDesigner />;
    case "sequences": return <FileManager />;
    case "pathways":
      return (
        <div className="flex items-center justify-center h-full text-gray-400">
          <div className="text-center space-y-2">
            <p className="text-2xl">🗺️</p>
            <p>Pathway Analysis</p>
            <p className="text-xs text-gray-600">Reactome, KEGG, STRING — search from the DB Search panel</p>
          </div>
        </div>
      );
    case "genome":
      return (
        <div className="flex items-center justify-center h-full text-gray-400">
          <div className="text-center space-y-2">
            <p className="text-2xl">🌐</p>
            <p>Genome Browser</p>
            <p className="text-xs text-gray-600">IGV.js integration — coming next</p>
          </div>
        </div>
      );
    default: return null;
  }
}

export default function App() {
  return (
    <QueryClientProvider client={queryClient}>
      <div className="flex h-screen bg-gray-950 text-gray-100 overflow-hidden">
        <Sidebar />
        <main className="flex-1 overflow-hidden">
          <MainPanel />
        </main>
      </div>
    </QueryClientProvider>
  );
}
