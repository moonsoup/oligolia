"""PCR primer design and restriction enzyme analysis panel."""

from __future__ import annotations
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, QPushButton,
    QLabel, QSpinBox, QDoubleSpinBox, QGroupBox, QTableWidget,
    QTableWidgetItem, QProgressBar, QTabWidget, QMessageBox, QHeaderView,
)
from PyQt6.QtGui import QColor

from backend.routers.primers import design_primers, restriction_sites, PrimerDesignRequest, RestrictionRequest
from ..workers import Worker


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

        btn_pcr = QPushButton("Design Primers")
        btn_pcr.setObjectName("primary")
        btn_pcr.clicked.connect(self._run_pcr)
        pcr_layout.addWidget(btn_pcr)

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

        layout.addWidget(tabs)

    def set_template(self, seq: str) -> None:
        self._template.setPlainText(seq)

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
