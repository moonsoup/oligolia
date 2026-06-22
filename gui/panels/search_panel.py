"""Database search panel — searches NCBI, Ensembl, UniProt, KEGG simultaneously."""

from __future__ import annotations
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton,
    QTableWidget, QTableWidgetItem, QCheckBox, QLabel, QGroupBox, QMessageBox, QHeaderView,
)
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtGui import QColor

from backend.services import NCBIClient, EnsemblClient, UniProtClient, KEGGClient
from backend.models.sequence import Sequence, MoleculeType
from ..workers import Worker


DB_COLORS = {
    "ncbi_gene": QColor("#1e3a5f"),
    "ensembl": QColor("#1e3a2e"),
    "uniprot": QColor("#3d2a00"),
    "kegg": QColor("#3d1a00"),
}


class SearchPanel(QWidget):
    sequence_fetched = pyqtSignal(object)  # emits Sequence

    def __init__(self) -> None:
        super().__init__()
        self._ncbi = NCBIClient()
        self._ensembl = EnsemblClient()
        self._uniprot = UniProtClient()
        self._kegg = KEGGClient()
        self._results: list[dict] = []
        self._worker: Worker | None = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        # Search bar
        search_row = QHBoxLayout()
        self._query = QLineEdit()
        self._query.setPlaceholderText("Search gene (e.g. BRCA2, TP53, CFTR)…")
        self._query.returnPressed.connect(self._run_search)
        search_row.addWidget(self._query, 1)
        self._species = QLineEdit("homo sapiens")
        self._species.setMaximumWidth(160)
        search_row.addWidget(QLabel("Species:"))
        search_row.addWidget(self._species)
        self._btn_search = QPushButton("Search")
        self._btn_search.setObjectName("primary")
        self._btn_search.clicked.connect(self._run_search)
        search_row.addWidget(self._btn_search)
        layout.addLayout(search_row)

        # Database checkboxes
        db_grp = QGroupBox("Databases")
        db_row = QHBoxLayout(db_grp)
        self._cb_ncbi = QCheckBox("NCBI Gene"); self._cb_ncbi.setChecked(True)
        self._cb_ensembl = QCheckBox("Ensembl"); self._cb_ensembl.setChecked(True)
        self._cb_uniprot = QCheckBox("UniProt"); self._cb_uniprot.setChecked(True)
        self._cb_kegg = QCheckBox("KEGG"); self._cb_kegg.setChecked(False)
        for cb in [self._cb_ncbi, self._cb_ensembl, self._cb_uniprot, self._cb_kegg]:
            db_row.addWidget(cb)
        db_row.addStretch()
        layout.addWidget(db_grp)

        # Status
        self._status = QLabel("Enter a gene name and click Search.")
        self._status.setObjectName("subheading")
        layout.addWidget(self._status)

        # Results table
        self._table = QTableWidget()
        self._table.setColumnCount(5)
        self._table.setHorizontalHeaderLabels(["Name", "Database", "Description", "Organism", "Accession"])
        self._table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self._table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.itemDoubleClicked.connect(self._on_double_click)
        self._table.currentCellChanged.connect(lambda row, *_: self._on_row_changed(row))
        layout.addWidget(self._table)

        # Detail / action row
        bottom = QHBoxLayout()
        self._detail = QLabel("Double-click a result to load its sequence into the editor.")
        self._detail.setObjectName("subheading")
        self._detail.setWordWrap(True)
        bottom.addWidget(self._detail, 1)
        self._btn_load = QPushButton("Load sequence")
        self._btn_load.setObjectName("primary")
        self._btn_load.setEnabled(False)
        self._btn_load.clicked.connect(self._load_selected)
        bottom.addWidget(self._btn_load)
        layout.addLayout(bottom)

    def _run_search(self) -> None:
        query = self._query.text().strip()
        if not query:
            return
        self._btn_search.setEnabled(False)
        self._btn_search.setText("Searching…")
        self._status.setText(f"Searching for '{query}'…")
        self._table.setRowCount(0)
        self._results.clear()

        dbs = []
        if self._cb_ncbi.isChecked():
            dbs.append("ncbi")
        if self._cb_ensembl.isChecked():
            dbs.append("ensembl")
        if self._cb_uniprot.isChecked():
            dbs.append("uniprot")
        if self._cb_kegg.isChecked():
            dbs.append("kegg")

        self._worker = Worker(self._do_search, query, self._species.text().strip(), dbs)
        self._worker.result.connect(self._on_search_done)
        self._worker.error.connect(self._on_search_error)
        self._worker.start()

    def _do_search(self, query: str, species: str, dbs: list[str]) -> dict:
        results = []
        errors: dict[str, str] = {}

        if "ncbi" in dbs:
            try:
                records = self._ncbi.search_genes(query, organism=species, max_results=10)
                for r in records:
                    results.append({
                        "name": r.get("name", ""),
                        "database": "ncbi_gene",
                        "description": r.get("description", ""),
                        "organism": r.get("organism", {}).get("scientificname", species),
                        "accession": str(r.get("uid", "")),
                        "url": f"https://www.ncbi.nlm.nih.gov/gene/{r.get('uid', '')}",
                        "_raw": r,
                    })
            except Exception as e:
                errors["NCBI"] = str(e)

        if "ensembl" in dbs:
            try:
                sp_slug = species.lower().replace(" ", "_")
                gene = self._ensembl.lookup_symbol(sp_slug, query)
                if gene and gene.get("id"):
                    results.append({
                        "name": gene.get("display_name", query),
                        "database": "ensembl",
                        "description": gene.get("description", ""),
                        "organism": gene.get("species", species),
                        "accession": gene.get("id", ""),
                        "url": f"https://www.ensembl.org/id/{gene.get('id', '')}",
                        "_raw": gene,
                    })
            except Exception:
                pass  # not found is normal

        if "uniprot" in dbs:
            try:
                data = self._uniprot.search_by_gene(query, organism=species, size=5)
                for entry in data.get("results", []):
                    acc = entry.get("primaryAccession", "")
                    genes = entry.get("genes", [{}])
                    gname = genes[0].get("geneName", {}).get("value", query) if genes else query
                    desc_obj = entry.get("proteinDescription", {}).get("recommendedName", {})
                    desc = desc_obj.get("fullName", {}).get("value", "") if desc_obj else ""
                    results.append({
                        "name": gname, "database": "uniprot",
                        "description": desc, "organism": species,
                        "accession": acc, "url": f"https://www.uniprot.org/uniprotkb/{acc}",
                        "_raw": entry,
                    })
            except Exception:
                pass

        if "kegg" in dbs:
            try:
                matches = self._kegg.search("genes", query)
                for kid, desc in matches[:5]:
                    results.append({
                        "name": kid, "database": "kegg",
                        "description": desc, "organism": species,
                        "accession": kid, "url": f"https://www.kegg.jp/entry/{kid}",
                        "_raw": {},
                    })
            except Exception:
                pass

        return {"results": results, "errors": errors}

    def _on_search_done(self, payload: dict) -> None:
        results = payload.get("results", [])
        errors = payload.get("errors", {})

        self._results = results
        self._table.setRowCount(len(results))
        for row, r in enumerate(results):
            for col, key in enumerate(["name", "database", "description", "organism", "accession"]):
                item = QTableWidgetItem(str(r.get(key, "")))
                bg = DB_COLORS.get(r.get("database", ""), QColor("#1e293b"))
                item.setBackground(bg)
                self._table.setItem(row, col, item)

        status = f"{len(results)} result{'s' if len(results) != 1 else ''} found."
        if errors:
            failed = ", ".join(f"{db} ({msg.split(chr(10))[0][:60]})" for db, msg in errors.items())
            status += f"  ⚠ Search failed for: {failed}"
        self._status.setText(status)

        self._btn_search.setEnabled(True)
        self._btn_search.setText("Search")

    def _on_search_error(self, err: str) -> None:
        self._status.setText(f"Error: {err}")
        self._btn_search.setEnabled(True)
        self._btn_search.setText("Search")

    def _on_row_changed(self, row: int) -> None:
        self._btn_load.setEnabled(0 <= row < len(self._results))

    def _on_double_click(self, item: QTableWidgetItem) -> None:
        self._load_selected()

    def _load_selected(self) -> None:
        row = self._table.currentRow()
        if row < 0 or row >= len(self._results):
            return
        result = self._results[row]
        self._btn_load.setEnabled(False)
        self._btn_load.setText("Fetching…")

        self._fetch_worker = Worker(self._fetch_sequence, result)
        self._fetch_worker.result.connect(self._on_fetch_done)
        self._fetch_worker.error.connect(lambda e: (
            QMessageBox.critical(self, "Fetch failed", e),
            self._btn_load.setEnabled(True),
            self._btn_load.setText("Load sequence"),
        ))
        self._fetch_worker.start()

    def _fetch_sequence(self, result: dict) -> Sequence:
        db = result.get("database", "")
        acc = result.get("accession", "")
        if db == "ncbi_gene" and acc:
            # acc is a Gene UID — must convert to nuccore via elink first
            fasta_text = self._ncbi.fetch_fasta_for_gene(acc)
            lines = fasta_text.splitlines()
            # FASTA header contains the actual accession (e.g. NM_000795.4)
            header = next((ln for ln in lines if ln.startswith(">")), "")
            nuccore_acc = header[1:].split()[0] if header else acc
            seq_str = "".join(ln for ln in lines if not ln.startswith(">"))
            return Sequence(
                id=nuccore_acc,
                name=result.get("name", nuccore_acc),
                description=result.get("description", header[1:] if header else ""),
                seq=seq_str.upper(),
                molecule_type=MoleculeType.DNA,
                accession=nuccore_acc,
                source_db="ncbi",
            )
        if db == "ensembl" and acc:
            data = self._ensembl.sequence_id(acc, seq_type="cdna")
            return Sequence(
                id=acc, name=result.get("name", acc),
                description=data.get("desc", ""),
                seq=data.get("seq", "").upper(), molecule_type=MoleculeType.DNA, accession=acc,
            )
        if db == "uniprot" and acc:
            fasta = self._uniprot.get_fasta(acc)
            seq_str = "".join(ln for ln in fasta.splitlines() if not ln.startswith(">"))
            return Sequence(
                id=acc, name=result.get("name", acc),
                description=result.get("description", ""),
                seq=seq_str.upper(), molecule_type=MoleculeType.PROTEIN, accession=acc,
            )
        raise ValueError(f"Cannot fetch sequence for database '{db}'")

    def _on_fetch_done(self, seq: Sequence) -> None:
        self.sequence_fetched.emit(seq)
        self._btn_load.setEnabled(True)
        self._btn_load.setText("Load sequence")
        self._status.setText(f"Loaded: {seq.id} ({len(seq.seq):,} bp)")
