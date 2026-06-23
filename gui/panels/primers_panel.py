"""PCR primer design and restriction enzyme analysis panel."""

from __future__ import annotations
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import json
from pathlib import Path

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, QPushButton,
    QLabel, QSpinBox, QDoubleSpinBox, QGroupBox, QTableWidget,
    QTableWidgetItem, QProgressBar, QTabWidget, QMessageBox, QHeaderView,
    QComboBox, QInputDialog,
)
from PyQt6.QtGui import QColor

from backend.routers.primers import design_primers, restriction_sites, PrimerDesignRequest, RestrictionRequest
from ..workers import Worker

# Presets stored in user config dir
_PRESETS_FILE = Path.home() / ".oligolia" / "primer_presets.json"
_BUILTIN_PRESETS = {
    "Standard PCR":      {"prod_min": 100,  "prod_max": 2000, "len_min": 18, "len_max": 22, "tm_min": 55.0, "tm_max": 65.0},
    "Long-range PCR":    {"prod_min": 2000, "prod_max": 15000,"len_min": 20, "len_max": 25, "tm_min": 60.0, "tm_max": 68.0},
    "Colony screening":  {"prod_min": 300,  "prod_max": 1500, "len_min": 18, "len_max": 22, "tm_min": 50.0, "tm_max": 60.0},
    "RT-qPCR":           {"prod_min": 80,   "prod_max": 200,  "len_min": 18, "len_max": 22, "tm_min": 58.0, "tm_max": 64.0},
    "Sequencing primer": {"prod_min": 100,  "prod_max": 1000, "len_min": 18, "len_max": 22, "tm_min": 55.0, "tm_max": 65.0},
}


class PrimersPanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self._worker: Worker | None = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        # Template
        tmpl_grp = QGroupBox("Template Sequence")
        tmpl_layout = QVBoxLayout(tmpl_grp)
        self._template = QTextEdit()
        self._template.setPlaceholderText("Paste or load template DNA sequence…")
        self._template.setMaximumHeight(70)
        tmpl_layout.addWidget(self._template)
        layout.addWidget(tmpl_grp)

        tabs = QTabWidget()

        # ── PCR Primer Design ────────────────────────────────────────────────
        pcr_widget = QWidget()
        pcr_layout = QVBoxLayout(pcr_widget)

        params_grp = QGroupBox("Parameters")
        params_layout = QHBoxLayout(params_grp)

        def spin(label: str, lo: int, hi: int, val: int, parent_layout=params_layout):
            parent_layout.addWidget(QLabel(label))
            s = QSpinBox(); s.setRange(lo, hi); s.setValue(val)
            parent_layout.addWidget(s)
            return s

        self._prod_min = spin("Product min (bp):", 50, 5000, 100)
        self._prod_max = spin("Product max (bp):", 50, 10000, 600)
        self._len_min = spin("Primer min:", 15, 35, 18)
        self._len_max = spin("Primer max:", 15, 35, 22)
        params_layout.addWidget(QLabel("Tm min:"))
        self._tm_min = QDoubleSpinBox(); self._tm_min.setRange(30, 90); self._tm_min.setValue(55.0)
        params_layout.addWidget(self._tm_min)
        params_layout.addWidget(QLabel("Tm max:"))
        self._tm_max = QDoubleSpinBox(); self._tm_max.setRange(30, 90); self._tm_max.setValue(65.0)
        params_layout.addWidget(self._tm_max)
        params_layout.addStretch()
        pcr_layout.addWidget(params_grp)

        # Preset row
        preset_row = QHBoxLayout()
        preset_row.addWidget(QLabel("Preset:"))
        self._preset_combo = QComboBox()
        self._preset_combo.setMinimumWidth(180)
        self._preset_combo.currentTextChanged.connect(self._load_preset)
        preset_row.addWidget(self._preset_combo)
        btn_save_preset = QPushButton("Save…")
        btn_save_preset.setObjectName("secondary")
        btn_save_preset.clicked.connect(self._save_preset)
        preset_row.addWidget(btn_save_preset)
        btn_del_preset = QPushButton("Delete")
        btn_del_preset.setObjectName("secondary")
        btn_del_preset.clicked.connect(self._delete_preset)
        preset_row.addWidget(btn_del_preset)
        preset_row.addStretch()

        btn_pcr = QPushButton("Design Primers")
        btn_pcr.setObjectName("primary")
        btn_pcr.clicked.connect(self._run_pcr)
        preset_row.addWidget(btn_pcr)
        pcr_layout.addLayout(preset_row)

        self._refresh_presets()

        self._pcr_progress = QProgressBar(); self._pcr_progress.setRange(0, 0); self._pcr_progress.hide()
        pcr_layout.addWidget(self._pcr_progress)
        self._pcr_status = QLabel(""); self._pcr_status.setObjectName("subheading")
        pcr_layout.addWidget(self._pcr_status)

        self._pcr_table = QTableWidget()
        self._pcr_table.setColumnCount(7)
        self._pcr_table.setHorizontalHeaderLabels(
            ["Pair", "Fwd sequence", "Rev sequence", "Product (bp)", "Fwd Tm", "Rev Tm", "Penalty"])
        self._pcr_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self._pcr_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self._pcr_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        pcr_layout.addWidget(self._pcr_table)
        tabs.addTab(pcr_widget, "PCR Primers")

        # ── Restriction Enzymes ───────────────────────────────────────────────
        re_widget = QWidget()
        re_layout = QVBoxLayout(re_widget)
        btn_re = QPushButton("Find Restriction Sites")
        btn_re.setObjectName("secondary")
        btn_re.clicked.connect(self._run_restriction)
        re_layout.addWidget(btn_re)
        self._re_status = QLabel(""); self._re_status.setObjectName("subheading")
        re_layout.addWidget(self._re_status)

        self._re_table = QTableWidget()
        self._re_table.setColumnCount(4)
        self._re_table.setHorizontalHeaderLabels(["Enzyme", "Recognition", "Count", "Positions"])
        self._re_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        self._re_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        re_layout.addWidget(self._re_table)
        tabs.addTab(re_widget, "Restriction Sites")

        # ── Restriction Digest ────────────────────────────────────────────────
        dig_widget = QWidget()
        dig_layout = QVBoxLayout(dig_widget)

        dig_ctrl = QHBoxLayout()
        dig_ctrl.addWidget(QLabel("Enzymes (comma-separated):"))
        self._dig_enzymes = QComboBox()
        self._dig_enzymes.setEditable(True)
        self._dig_enzymes.setMinimumWidth(260)
        from backend.routers.primers import RESTRICTION_ENZYMES as _RE
        common = ["EcoRI,BamHI", "EcoRI,HindIII", "BamHI,HindIII", "EcoRI,BamHI,HindIII",
                  "NcoI,XhoI", "NdeI,XhoI", "NdeI,BamHI", "SalI,XbaI"]
        self._dig_enzymes.addItems(common)
        self._dig_enzymes.setCurrentText("")
        dig_ctrl.addWidget(self._dig_enzymes, 1)
        btn_dig = QPushButton("Simulate Digest")
        btn_dig.setObjectName("primary")
        btn_dig.clicked.connect(self._run_digest)
        dig_ctrl.addWidget(btn_dig)
        dig_layout.addLayout(dig_ctrl)

        self._dig_status = QLabel(""); self._dig_status.setObjectName("subheading")
        dig_layout.addWidget(self._dig_status)

        self._dig_table = QTableWidget()
        self._dig_table.setColumnCount(4)
        self._dig_table.setHorizontalHeaderLabels(["Fragment", "Size (bp)", "Start", "End"])
        self._dig_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._dig_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        dig_layout.addWidget(self._dig_table)
        tabs.addTab(dig_widget, "Digest")

        layout.addWidget(tabs)

    def set_template(self, seq: str) -> None:
        self._template.setPlainText(seq)

    # ── Presets ───────────────────────────────────────────────────────────────

    def _all_presets(self) -> dict:
        user = {}
        if _PRESETS_FILE.exists():
            try:
                user = json.loads(_PRESETS_FILE.read_text())
            except Exception:
                pass
        return {**_BUILTIN_PRESETS, **user}

    def _refresh_presets(self) -> None:
        self._preset_combo.blockSignals(True)
        current = self._preset_combo.currentText()
        self._preset_combo.clear()
        self._preset_combo.addItem("— select preset —")
        for name in self._all_presets():
            self._preset_combo.addItem(name)
        idx = self._preset_combo.findText(current)
        self._preset_combo.setCurrentIndex(max(0, idx))
        self._preset_combo.blockSignals(False)

    def _load_preset(self, name: str) -> None:
        p = self._all_presets().get(name)
        if not p:
            return
        self._prod_min.setValue(p["prod_min"])
        self._prod_max.setValue(p["prod_max"])
        self._len_min.setValue(p["len_min"])
        self._len_max.setValue(p["len_max"])
        self._tm_min.setValue(p["tm_min"])
        self._tm_max.setValue(p["tm_max"])

    def _save_preset(self) -> None:
        name, ok = QInputDialog.getText(self, "Save Preset", "Preset name:")
        if not ok or not name.strip():
            return
        name = name.strip()
        if name in _BUILTIN_PRESETS:
            QMessageBox.warning(self, "Reserved name", f"'{name}' is a built-in preset — choose a different name.")
            return
        _PRESETS_FILE.parent.mkdir(parents=True, exist_ok=True)
        user = {}
        if _PRESETS_FILE.exists():
            try:
                user = json.loads(_PRESETS_FILE.read_text())
            except Exception:
                pass
        user[name] = {
            "prod_min": self._prod_min.value(), "prod_max": self._prod_max.value(),
            "len_min": self._len_min.value(),   "len_max": self._len_max.value(),
            "tm_min": self._tm_min.value(),      "tm_max": self._tm_max.value(),
        }
        _PRESETS_FILE.write_text(json.dumps(user, indent=2))
        self._refresh_presets()
        self._preset_combo.setCurrentText(name)

    def _delete_preset(self) -> None:
        name = self._preset_combo.currentText()
        if name in _BUILTIN_PRESETS or name == "— select preset —":
            QMessageBox.warning(self, "Cannot delete", "Built-in presets cannot be deleted.")
            return
        if not _PRESETS_FILE.exists():
            return
        try:
            user = json.loads(_PRESETS_FILE.read_text())
            user.pop(name, None)
            _PRESETS_FILE.write_text(json.dumps(user, indent=2))
        except Exception:
            pass
        self._refresh_presets()

    def _run_pcr(self) -> None:
        template = self._template.toPlainText().strip().upper().replace(" ", "").replace("\n", "")
        if not template:
            QMessageBox.warning(self, "No template", "Paste a template sequence first.")
            return
        self._pcr_progress.show()
        self._pcr_table.setRowCount(0)
        req = PrimerDesignRequest(
            template=template,
            product_min=self._prod_min.value(),
            product_max=self._prod_max.value(),
            primer_len_min=self._len_min.value(),
            primer_len_max=self._len_max.value(),
            tm_min=self._tm_min.value(),
            tm_max=self._tm_max.value(),
            max_pairs=10,
        )
        self._worker = Worker(design_primers, req)
        self._worker.result.connect(self._on_pcr_done)
        self._worker.error.connect(lambda e: (self._pcr_progress.hide(), self._pcr_status.setText(f"Error: {e}")))
        self._worker.start()

    def _on_pcr_done(self, pairs: list) -> None:
        self._pcr_progress.hide()
        self._pcr_status.setText(f"{len(pairs)} primer pair{'s' if len(pairs) != 1 else ''} found.")
        self._pcr_table.setRowCount(len(pairs))
        for i, pair in enumerate(pairs):
            fwd, rev = pair.forward, pair.reverse
            vals = [
                str(i + 1), fwd.sequence, rev.sequence,
                str(pair.product_size),
                f"{fwd.tm:.1f}°C", f"{rev.tm:.1f}°C",
                f"{pair.penalty:.3f}",
            ]
            for col, val in enumerate(vals):
                item = QTableWidgetItem(val)
                if col in (1, 2):
                    from PyQt6.QtGui import QFont
                    item.setFont(QFont("JetBrains Mono,Courier New", 11))
                self._pcr_table.setItem(i, col, item)

    def _run_digest(self) -> None:
        template = self._template.toPlainText().strip().upper().replace(" ", "").replace("\n", "")
        if not template:
            QMessageBox.warning(self, "No template", "Paste a template sequence first.")
            return
        enzyme_text = self._dig_enzymes.currentText().strip()
        if not enzyme_text:
            QMessageBox.warning(self, "No enzymes", "Enter one or more enzyme names.")
            return
        enzymes = [e.strip() for e in enzyme_text.replace("+", ",").split(",") if e.strip()]
        from backend.routers.primers import digest, DigestRequest
        try:
            result = digest(DigestRequest(template=template, enzymes=enzymes))
            n = len(result.fragments)
            cuts = len(result.cut_positions)
            self._dig_status.setText(
                f"{cuts} cut site{'s' if cuts != 1 else ''} → {n} fragment{'s' if n != 1 else ''}  "
                f"(template: {result.template_length:,} bp)"
            )
            self._dig_table.setRowCount(n)
            for i, frag in enumerate(result.fragments):
                vals = [str(i + 1), f"{frag.length:,}", str(frag.start), str(frag.end)]
                for col, val in enumerate(vals):
                    item = QTableWidgetItem(val)
                    if i % 2 == 0:
                        item.setBackground(QColor("#1a2030"))
                    self._dig_table.setItem(i, col, item)
        except Exception as e:
            self._dig_status.setText(f"Error: {e}")

    def _run_restriction(self) -> None:
        template = self._template.toPlainText().strip().upper().replace(" ", "").replace("\n", "")
        if not template:
            return
        req = RestrictionRequest(template=template)
        try:
            sites = restriction_sites(req)
            self._re_status.setText(f"{len(sites)} enzyme{'s' if len(sites) != 1 else ''} cut.")
            self._re_table.setRowCount(len(sites))
            for i, s in enumerate(sites):
                vals = [s.enzyme, s.cut_pattern, str(s.count), ", ".join(str(p) for p in s.positions[:10])]
                for col, val in enumerate(vals):
                    item = QTableWidgetItem(val)
                    if s.count >= 2:
                        item.setBackground(QColor("#1e3a2e"))
                    self._re_table.setItem(i, col, item)
        except Exception as e:
            self._re_status.setText(f"Error: {e}")
