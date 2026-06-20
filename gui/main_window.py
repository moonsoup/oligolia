"""Oligolia main application window."""

from __future__ import annotations
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QTabWidget,
    QMessageBox,
    QFileDialog, QApplication,
)
from PyQt6.QtGui import QAction

from .styles import DARK_STYLESHEET
from .panels import (
    SequencePanel, SearchPanel, CRISPRPanel,
    AlignmentPanel, PrimersPanel, VariantsPanel,
)
from backend.formats import write_fasta, write_genbank
from backend.models.sequence import Sequence


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Oligolia — Gene Editing Platform")
        self.resize(1400, 900)
        self.setMinimumSize(900, 600)
        self._build_menu()
        self._build_ui()
        self._build_status()
        self.setStyleSheet(DARK_STYLESHEET)

    # ── Menu ─────────────────────────────────────────────────────────────────

    def _build_menu(self) -> None:
        mb = self.menuBar()

        # File
        file_menu = mb.addMenu("File")
        act_open = QAction("Open sequence file…", self)
        act_open.setShortcut("Ctrl+O")
        act_open.triggered.connect(self._open_file_menu)
        file_menu.addAction(act_open)

        act_open_vcf = QAction("Open VCF file…", self)
        act_open_vcf.triggered.connect(lambda: self._tabs.setCurrentWidget(self._variants_panel) or self._variants_panel._open_vcf())
        file_menu.addAction(act_open_vcf)

        file_menu.addSeparator()
        act_save_fasta = QAction("Save all sequences as FASTA…", self)
        act_save_fasta.setShortcut("Ctrl+S")
        act_save_fasta.triggered.connect(self._save_fasta)
        file_menu.addAction(act_save_fasta)

        act_save_gb = QAction("Save all sequences as GenBank…", self)
        act_save_gb.triggered.connect(self._save_genbank)
        file_menu.addAction(act_save_gb)

        file_menu.addSeparator()
        act_quit = QAction("Quit", self)
        act_quit.setShortcut("Ctrl+Q")
        act_quit.triggered.connect(QApplication.quit)
        file_menu.addAction(act_quit)

        # Tools
        tools_menu = mb.addMenu("Tools")
        act_gc = QAction("GC Content of active sequence", self)
        act_gc.triggered.connect(self._show_gc)
        tools_menu.addAction(act_gc)

        act_codon = QAction("Codon Usage of active sequence", self)
        act_codon.triggered.connect(self._show_codon_usage)
        tools_menu.addAction(act_codon)

        # Help
        help_menu = mb.addMenu("Help")
        act_about = QAction("About Oligolia", self)
        act_about.triggered.connect(self._show_about)
        help_menu.addAction(act_about)

    # ── UI ───────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        central = QWidget()
        self.setCentralWidget(central)
        layout = QHBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self._tabs = QTabWidget()
        self._tabs.setTabPosition(QTabWidget.TabPosition.West)
        self._tabs.setMovable(True)

        # Instantiate panels
        self._seq_panel = SequencePanel()
        self._search_panel = SearchPanel()
        self._crispr_panel = CRISPRPanel()
        self._align_panel = AlignmentPanel()
        self._primers_panel = PrimersPanel()
        self._variants_panel = VariantsPanel()

        # Connect cross-panel signals
        self._seq_panel.sequence_selected.connect(self._on_seq_selected)
        self._search_panel.sequence_fetched.connect(self._seq_panel.add_sequence)
        self._search_panel.sequence_fetched.connect(
            lambda seq: (self._tabs.setCurrentWidget(self._seq_panel), self._status.showMessage(f"Loaded: {seq.id}"))
        )

        def _add(panel: QWidget, icon: str, label: str) -> None:
            self._tabs.addTab(panel, f"{icon} {label}")

        _add(self._seq_panel,     "🧬", "Sequences")
        _add(self._search_panel,  "🔍", "DB Search")
        _add(self._crispr_panel,  "✂️", "CRISPR")
        _add(self._align_panel,   "↔", "Alignment")
        _add(self._primers_panel, "🔩", "Primers")
        _add(self._variants_panel,"🔬", "Variants")

        layout.addWidget(self._tabs)

    def _build_status(self) -> None:
        self._status = self.statusBar()
        self._status.showMessage("Oligolia ready — open a file or search a database to start.")

    # ── Signal handlers ───────────────────────────────────────────────────────

    def _on_seq_selected(self, seq: Sequence) -> None:
        self._status.showMessage(
            f"Active: {seq.id}  |  {seq.molecule_type.value}  |  {seq.length:,} bp"
        )
        # Pre-fill other panels with the active sequence
        if seq.molecule_type.value in ("DNA", "RNA"):
            self._crispr_panel.set_target(seq.seq)
            self._primers_panel.set_template(seq.seq)

    # ── Menu actions ──────────────────────────────────────────────────────────

    def _open_file_menu(self) -> None:
        self._tabs.setCurrentWidget(self._seq_panel)
        self._seq_panel._open_file()

    def _save_fasta(self) -> None:
        seqs = list(self._seq_panel.sequences.values())
        if not seqs:
            QMessageBox.information(self, "Nothing to save", "Load or create sequences first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save FASTA", "", "FASTA (*.fa *.fasta);;All files (*)")
        if path:
            with open(path, "w") as f:
                f.write(write_fasta(seqs))
            self._status.showMessage(f"Saved {len(seqs)} sequences to {path}")

    def _save_genbank(self) -> None:
        seqs = list(self._seq_panel.sequences.values())
        if not seqs:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save GenBank", "", "GenBank (*.gb);;All files (*)")
        if path:
            with open(path, "w") as f:
                f.write(write_genbank(seqs))
            self._status.showMessage(f"Saved {len(seqs)} sequences to {path}")

    def _show_gc(self) -> None:
        if not self._seq_panel._active:
            QMessageBox.information(self, "No sequence", "Select a sequence first.")
            return
        seq = self._seq_panel._active.seq.upper()
        gc = (seq.count("G") + seq.count("C")) / len(seq) * 100 if seq else 0.0
        QMessageBox.information(
            self, "GC Content",
            f"Sequence: {self._seq_panel._active.id}\n"
            f"Length: {len(seq):,} bp\n"
            f"GC content: {gc:.2f}%\n"
            f"AT content: {100 - gc:.2f}%"
        )

    def _show_codon_usage(self) -> None:
        if not self._seq_panel._active:
            return
        seq = self._seq_panel._active.seq.upper()
        counts: dict[str, int] = {}
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if len(codon) == 3:
                counts[codon] = counts.get(codon, 0) + 1
        total = sum(counts.values())
        lines = [f"Codon usage for {self._seq_panel._active.id} ({total} codons total)\n"]
        for codon, cnt in sorted(counts.items()):
            pct = cnt / total * 100 if total else 0
            lines.append(f"{codon}  {cnt:4d}  {pct:.1f}%")
        msg = QMessageBox(self)
        msg.setWindowTitle("Codon Usage")
        msg.setText("\n".join(lines))
        msg.exec()

    def _show_about(self) -> None:
        QMessageBox.about(
            self, "About Oligolia",
            "<h2>Oligolia</h2>"
            "<p>Advanced gene editing and viewing platform.</p>"
            "<p><b>Features:</b><br>"
            "• Sequence editing (insert/delete/replace/RC/translate/transcribe)<br>"
            "• Multi-database search (NCBI, Ensembl, UniProt, KEGG)<br>"
            "• CRISPR guide RNA design (SpCas9, Cas12a, Cas13)<br>"
            "• Pairwise &amp; multiple sequence alignment<br>"
            "• PCR primer design + restriction enzyme analysis<br>"
            "• VCF variant viewer<br>"
            "• FASTA, FASTQ, GenBank, EMBL, GFF3/GTF, VCF support</p>"
            "<p><b>Version:</b> 0.1.0<br>"
            "<b>Backend:</b> Biopython 1.85</p>"
        )
