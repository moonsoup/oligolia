"""Pathway analysis panel — Reactome enrichment, KEGG, STRING interactions."""

from __future__ import annotations
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QTableWidget, QTableWidgetItem, QPlainTextEdit, QTabWidget,
    QHeaderView, QMessageBox, QComboBox, QSpinBox, QGroupBox,
)
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtGui import QColor, QDesktopServices
from PyQt6.QtCore import QUrl

from backend.services import ReactomeClient, KEGGClient, STRINGClient
from ..workers import Worker


class PathwaysPanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self._reactome = ReactomeClient()
        self._kegg = KEGGClient()
        self._string = STRINGClient()
        self._worker: Worker | None = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        # ── Input ─────────────────────────────────────────────────────────────
        input_grp = QGroupBox("Gene / Protein Input")
        input_layout = QVBoxLayout(input_grp)

        hint = QLabel("Enter gene symbols or UniProt IDs, one per line (e.g. BRCA1, TP53, EGFR)")
        hint.setObjectName("subheading")
        input_layout.addWidget(hint)

        self._gene_input = QPlainTextEdit()
        self._gene_input.setPlaceholderText("BRCA1\nTP53\nEGFR\nPTEN\nMYC")
        self._gene_input.setMaximumHeight(110)
        self._gene_input.setFont(self._gene_input.font())
        input_layout.addWidget(self._gene_input)

        ctrl_row = QHBoxLayout()
        ctrl_row.addWidget(QLabel("Species:"))
        self._species = QComboBox()
        self._species.addItems(["Homo sapiens", "Mus musculus", "Rattus norvegicus"])
        ctrl_row.addWidget(self._species)
        ctrl_row.addStretch()

        self._btn_reactome = QPushButton("Reactome Enrichment")
        self._btn_reactome.setObjectName("primary")
        self._btn_reactome.clicked.connect(self._run_reactome)
        ctrl_row.addWidget(self._btn_reactome)

        self._btn_string = QPushButton("STRING Network")
        self._btn_string.setObjectName("secondary")
        self._btn_string.clicked.connect(self._run_string)
        ctrl_row.addWidget(self._btn_string)

        input_layout.addLayout(ctrl_row)
        layout.addWidget(input_grp)

        self._status = QLabel("")
        self._status.setObjectName("subheading")
        layout.addWidget(self._status)

        # ── Results tabs ──────────────────────────────────────────────────────
        tabs = QTabWidget()

        # Reactome
        reactome_widget = QWidget()
        rl = QVBoxLayout(reactome_widget)
        self._reactome_table = QTableWidget()
        self._reactome_table.setColumnCount(5)
        self._reactome_table.setHorizontalHeaderLabels(
            ["Pathway", "Entities Found", "Entities Total", "p-value", "FDR"]
        )
        self._reactome_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self._reactome_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._reactome_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self._reactome_table.itemDoubleClicked.connect(self._open_reactome_pathway)
        rl.addWidget(self._reactome_table)
        rl.addWidget(QLabel("Double-click a pathway to open it in Reactome"))
        tabs.addTab(reactome_widget, "Reactome Enrichment")

        # STRING
        string_widget = QWidget()
        sl = QVBoxLayout(string_widget)
        self._string_table = QTableWidget()
        self._string_table.setColumnCount(4)
        self._string_table.setHorizontalHeaderLabels(
            ["Protein A", "Protein B", "Score", "Interaction Type"]
        )
        self._string_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._string_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._string_url: str = ""
        sl.addWidget(self._string_table)
        self._btn_string_open = QPushButton("Open full network in STRING")
        self._btn_string_open.setEnabled(False)
        self._btn_string_open.clicked.connect(self._open_string_network)
        sl.addWidget(self._btn_string_open)
        tabs.addTab(string_widget, "STRING Network")

        layout.addWidget(tabs)

    def _gene_list(self) -> list[str]:
        return [g.strip() for g in self._gene_input.toPlainText().splitlines() if g.strip()]

    # ── Reactome ──────────────────────────────────────────────────────────────

    def _run_reactome(self) -> None:
        genes = self._gene_list()
        if not genes:
            QMessageBox.warning(self, "No input", "Enter at least one gene symbol.")
            return
        self._btn_reactome.setEnabled(False)
        self._status.setText(f"Running Reactome enrichment for {len(genes)} genes…")
        self._reactome_table.setRowCount(0)
        species = self._species.currentText()
        self._worker = Worker(self._fetch_reactome, genes, species)
        self._worker.result.connect(self._on_reactome_done)
        self._worker.error.connect(lambda e: (
            self._status.setText(f"Reactome error: {e}"),
            self._btn_reactome.setEnabled(True),
        ))
        self._worker.start()

    def _fetch_reactome(self, genes: list[str], species: str) -> list[dict]:
        result = self._reactome.mapping_analysis(genes, species=species)
        pathways = result.get("pathways", result.get("entities", []))
        if isinstance(pathways, dict):
            pathways = pathways.get("pathways", [])
        return pathways if isinstance(pathways, list) else []

    def _on_reactome_done(self, pathways: list[dict]) -> None:
        self._btn_reactome.setEnabled(True)
        self._reactome_table.setRowCount(len(pathways))
        for i, p in enumerate(pathways):
            name = p.get("name", p.get("pathway", {}).get("name", ""))
            found = str(p.get("entities", {}).get("found", p.get("found", "")))
            total = str(p.get("entities", {}).get("total", p.get("total", "")))
            pval = f"{p.get('entities', {}).get('pValue', p.get('pValue', '')):.2e}" if p.get('entities', p).get('pValue') else ""
            fdr = f"{p.get('entities', {}).get('fdr', p.get('fdr', '')):.2e}" if p.get('entities', p).get('fdr') else ""
            for col, val in enumerate([name, found, total, pval, fdr]):
                item = QTableWidgetItem(str(val))
                if i % 2 == 0:
                    item.setBackground(QColor("#1a2030"))
                self._reactome_table.setItem(i, col, item)
            self._reactome_table.item(i, 0).setData(256, p.get("stId", p.get("pathway", {}).get("stId", "")))
        self._status.setText(f"{len(pathways)} enriched pathways found.")

    def _open_reactome_pathway(self, item: QTableWidgetItem) -> None:
        st_id = self._reactome_table.item(item.row(), 0).data(256)
        if st_id:
            QDesktopServices.openUrl(QUrl(f"https://reactome.org/PathwayBrowser/#/{st_id}"))

    # ── STRING ────────────────────────────────────────────────────────────────

    def _run_string(self) -> None:
        genes = self._gene_list()
        if not genes:
            QMessageBox.warning(self, "No input", "Enter at least one gene symbol.")
            return
        self._btn_string.setEnabled(False)
        self._status.setText(f"Fetching STRING interactions for {len(genes)} proteins…")
        self._string_table.setRowCount(0)
        self._worker = Worker(self._fetch_string, genes)
        self._worker.result.connect(self._on_string_done)
        self._worker.error.connect(lambda e: (
            self._status.setText(f"STRING error: {e}"),
            self._btn_string.setEnabled(True),
        ))
        self._worker.start()

    def _fetch_string(self, genes: list[str]) -> list[dict]:
        return self._string.network(genes, species=9606, required_score=400)

    def _on_string_done(self, interactions: list[dict]) -> None:
        self._btn_string.setEnabled(True)
        self._string_table.setRowCount(len(interactions))
        for i, ix in enumerate(interactions):
            a = ix.get("preferredName_A", ix.get("stringId_A", ""))
            b = ix.get("preferredName_B", ix.get("stringId_B", ""))
            score = f"{ix.get('score', ix.get('combined_score', 0)):.3f}"
            itype = ", ".join(k for k, v in ix.items()
                              if k in ("coexpression", "cooccurence", "fusion",
                                       "neighborhood", "experimental", "textmining")
                              and v and float(v) > 0.3)
            for col, val in enumerate([a, b, score, itype]):
                self._string_table.setItem(i, col, QTableWidgetItem(str(val)))

        genes = self._gene_list()
        self._string_url = (
            f"https://string-db.org/network/{','.join(genes)}"
        )
        self._btn_string_open.setEnabled(True)
        self._status.setText(f"{len(interactions)} interactions found.")

    def _open_string_network(self) -> None:
        if self._string_url:
            QDesktopServices.openUrl(QUrl(self._string_url))

    def load_genes(self, genes: list[str]) -> None:
        """Pre-populate the gene list (called from other panels)."""
        self._gene_input.setPlainText("\n".join(genes))
