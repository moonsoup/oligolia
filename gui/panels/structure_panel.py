"""Protein structure lookup/prediction + putative interaction points.

Data-only in the core app: no in-app 3D rendering (that lives in the separate
"Oligolia Structure Viewer" companion app to keep the core installer free of a
bundled Chromium/WebEngine runtime — see gui/plugins/structure_viewer_launcher.py).
This panel shows a source badge, a heuristic interaction-points table, and lets
the user save the .pdb, hand off to the companion viewer if installed, or flag a
prediction as wrong.
"""

from __future__ import annotations
import sys
import os
import json
import tempfile
import webbrowser
from urllib.parse import quote

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton,
    QTableWidget, QTableWidgetItem, QLabel, QGroupBox, QMessageBox,
    QHeaderView, QFileDialog,
)
from PyQt6.QtGui import QColor

from backend.models.structure import StructureRequest, StructureResult, InteractionPointsRequest
from backend.routers.structure import get_or_predict_structure, interaction_points as compute_interaction_points_endpoint
from ..workers import Worker
from ..plugins.structure_viewer_launcher import find_structure_viewer

try:
    from version import GITHUB_OWNER, GITHUB_REPO
except ImportError:
    GITHUB_OWNER, GITHUB_REPO = "moonsoup", "oligolia"

SOURCE_LABELS = {
    "experimental_pdb": "Experimental structure (PDB)",
    "predicted_alphafold_db": "Predicted structure (AlphaFold DB)",
    "predicted_esmfold": "Predicted structure (ESMFold, novel sequence)",
}

CLASS_COLORS = {
    "acidic": QColor("#4a1e1e"),
    "basic": QColor("#1e2a4a"),
    "polar": QColor("#1e3a2e"),
    "hydrophobic": QColor("#1e293b"),
}

DISCLAIMER = (
    "Putative interaction points are a heuristic surface-exposure + charge/polarity "
    "annotation, NOT a validated binding-site or docking prediction."
)


class StructurePanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self._sequence: str = ""
        self._gene_symbol: str | None = None
        self._uniprot_id: str | None = None
        self._result: StructureResult | None = None
        self._points: list[dict] = []
        self._worker: Worker | None = None
        self._points_worker: Worker | None = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        target_grp = QGroupBox("Target Protein")
        target_layout = QHBoxLayout(target_grp)
        self._target_label = QLabel("No protein sequence loaded — select one in Sequences.")
        self._target_label.setWordWrap(True)
        target_layout.addWidget(self._target_label, 1)
        self._gene_input = QLineEdit()
        self._gene_input.setPlaceholderText("Gene symbol (optional)")
        self._gene_input.setMaximumWidth(160)
        target_layout.addWidget(self._gene_input)
        self._uniprot_input = QLineEdit()
        self._uniprot_input.setPlaceholderText("UniProt accession (optional)")
        self._uniprot_input.setMaximumWidth(180)
        target_layout.addWidget(self._uniprot_input)
        layout.addWidget(target_grp)

        action_row = QHBoxLayout()
        self._btn_predict = QPushButton("Get / Predict Structure")
        self._btn_predict.setObjectName("primary")
        self._btn_predict.setEnabled(False)
        self._btn_predict.clicked.connect(self._run_predict)
        action_row.addWidget(self._btn_predict)
        action_row.addStretch()
        layout.addLayout(action_row)

        self._badge = QLabel("")
        self._badge.setObjectName("subheading")
        layout.addWidget(self._badge)

        self._table = QTableWidget()
        self._table.setColumnCount(5)
        self._table.setHorizontalHeaderLabels(
            ["Residue #", "Residue", "Chain", "Classification", "Relative SASA"]
        )
        self._table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        self._table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSortingEnabled(True)
        layout.addWidget(self._table)

        disclaimer = QLabel(DISCLAIMER)
        disclaimer.setWordWrap(True)
        disclaimer.setObjectName("subheading")
        layout.addWidget(disclaimer)

        bottom = QHBoxLayout()
        self._btn_save = QPushButton("Save Structure (.pdb)")
        self._btn_save.setEnabled(False)
        self._btn_save.clicked.connect(self._save_pdb)
        bottom.addWidget(self._btn_save)

        self._btn_viewer = QPushButton("Open in 3D Viewer")
        self._btn_viewer.setEnabled(False)
        self._btn_viewer.clicked.connect(self._open_in_viewer)
        bottom.addWidget(self._btn_viewer)

        self._btn_flag = QPushButton("Flag as wrong")
        self._btn_flag.setEnabled(False)
        self._btn_flag.clicked.connect(self._flag_as_wrong)
        bottom.addWidget(self._btn_flag)

        bottom.addStretch()
        layout.addLayout(bottom)

    # ── Cross-panel wiring ───────────────────────────────────────────────────

    def set_target(self, sequence: str, gene_symbol: str | None = None, uniprot_id: str | None = None) -> None:
        self._sequence = sequence
        self._gene_symbol = gene_symbol
        self._uniprot_id = uniprot_id
        if uniprot_id:
            self._uniprot_input.setText(uniprot_id)
        if gene_symbol:
            self._gene_input.setText(gene_symbol)
        self._target_label.setText(f"Loaded protein: {len(sequence):,} residues")
        self._btn_predict.setEnabled(bool(sequence))

    # ── Predict / lookup ─────────────────────────────────────────────────────

    def _run_predict(self) -> None:
        self._btn_predict.setEnabled(False)
        self._btn_predict.setText("Working…")
        self._badge.setText("Looking up / predicting structure…")
        self._table.setRowCount(0)
        self._btn_save.setEnabled(False)
        self._btn_viewer.setEnabled(False)
        self._btn_flag.setEnabled(False)

        req = StructureRequest(
            sequence=self._sequence,
            gene_symbol=self._gene_input.text().strip() or None,
            uniprot_id=self._uniprot_input.text().strip() or None,
        )
        self._worker = Worker(get_or_predict_structure, req)
        self._worker.result.connect(self._on_predict_done)
        self._worker.error.connect(self._on_predict_error)
        self._worker.start()

    def _on_predict_done(self, result: StructureResult) -> None:
        self._result = result
        label = SOURCE_LABELS.get(result.source.value, result.source.value)
        badge_text = label
        if result.pdb_id:
            badge_text += f" — {result.pdb_id}"
        if result.confidence_note:
            badge_text += f"\n{result.confidence_note}"
        if result.warnings:
            badge_text += "\n" + "\n".join(f"⚠ {w}" for w in result.warnings)
        self._badge.setText(badge_text)

        self._btn_predict.setEnabled(True)
        self._btn_predict.setText("Get / Predict Structure")
        self._btn_save.setEnabled(True)
        self._btn_viewer.setEnabled(True)
        self._btn_flag.setEnabled(True)

        self._points_worker = Worker(
            compute_interaction_points_endpoint, InteractionPointsRequest(pdb_text=result.pdb_text)
        )
        self._points_worker.result.connect(self._on_points_done)
        self._points_worker.error.connect(lambda e: self._badge.setText(self._badge.text() + f"\n(interaction points failed: {e})"))
        self._points_worker.start()

    def _on_predict_error(self, err: str) -> None:
        self._btn_predict.setEnabled(True)
        self._btn_predict.setText("Get / Predict Structure")
        if "400" in err or "residues" in err.lower():
            self._badge.setText(f"⚠ {err}")
        elif any(term in err.lower() for term in ("connection", "network", "resolve", "timeout")):
            self._badge.setText("⚠ Could not reach the structure prediction service — check your connection.")
        else:
            self._badge.setText(f"⚠ No structure found and prediction failed: {err}")

    def _on_points_done(self, result) -> None:
        self._points = [p.model_dump() for p in result.points]
        self._table.setSortingEnabled(False)
        self._table.setRowCount(len(self._points))
        for row, p in enumerate(self._points):
            values = [
                str(p["residue_index"]), p["residue_name"], p["chain"],
                p["classification"] + (" ★" if p["is_putative_interaction_point"] else ""),
                f"{p['relative_sasa']:.2f}",
            ]
            for col, val in enumerate(values):
                item = QTableWidgetItem(val)
                if p["is_putative_interaction_point"]:
                    item.setBackground(CLASS_COLORS.get(p["classification"], QColor("#1e293b")))
                self._table.setItem(row, col, item)
        self._table.setSortingEnabled(True)

    # ── Actions ──────────────────────────────────────────────────────────────

    def _save_pdb(self) -> None:
        if not self._result:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save Structure", "", "PDB (*.pdb);;All files (*)")
        if path:
            with open(path, "w") as f:
                f.write(self._result.pdb_text)

    def _open_in_viewer(self) -> None:
        if not self._result:
            return
        viewer = find_structure_viewer()
        if viewer is None:
            box = QMessageBox(self)
            box.setWindowTitle("Structure Viewer not installed")
            box.setText(
                "The optional 3D structure viewer isn't installed.\n\n"
                "You can download it separately, or save the .pdb file and open it in "
                "an external viewer (PyMOL, ChimeraX, etc.)."
            )
            open_btn = box.addButton("Download Structure Viewer", QMessageBox.ButtonRole.ActionRole)
            box.addButton(QMessageBox.StandardButton.Close)
            box.exec()
            if box.clickedButton() == open_btn:
                webbrowser.open(f"https://github.com/{GITHUB_OWNER}/{GITHUB_REPO}/releases")
            return

        import subprocess
        tmp_dir = tempfile.mkdtemp(prefix="oligolia_structure_")
        pdb_path = os.path.join(tmp_dir, "structure.pdb")
        with open(pdb_path, "w") as f:
            f.write(self._result.pdb_text)
        args = [str(viewer), "--pdb", pdb_path]
        if self._points:
            points_path = os.path.join(tmp_dir, "interaction_points.json")
            with open(points_path, "w") as f:
                json.dump(self._points, f)
            args += ["--interactions", points_path]
        subprocess.Popen(args)

    def _flag_as_wrong(self) -> None:
        if not self._result:
            return
        accession = self._uniprot_input.text().strip() or self._gene_input.text().strip() or "unknown"
        source = SOURCE_LABELS.get(self._result.source.value, self._result.source.value)
        title = f"Structure prediction quality: {accession}"
        body = (
            f"**Source:** {source}\n"
            f"**Accession/gene:** {accession}\n"
            f"**PDB ID (if experimental):** {self._result.pdb_id or 'n/a'}\n"
            f"**Confidence note:** {self._result.confidence_note}\n\n"
            "**What looks wrong?**\n(describe here)\n"
        )
        url = (
            f"https://github.com/{GITHUB_OWNER}/{GITHUB_REPO}/issues/new"
            f"?title={quote(title)}&labels=structure-feedback&body={quote(body)}"
        )
        webbrowser.open(url)
