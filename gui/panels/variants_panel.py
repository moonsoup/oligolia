"""VCF variant viewer panel."""

from __future__ import annotations
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QTableWidget, QTableWidgetItem, QLineEdit, QFileDialog,
    QGroupBox, QTextEdit, QHeaderView, QMessageBox,
)
from PyQt6.QtGui import QColor

from backend.formats import parse_vcf
from backend.models.variant import Variant, VariantType


TYPE_COLORS = {
    VariantType.SNP: QColor("#1a2d3a"),
    VariantType.DEL: QColor("#3d1a1a"),
    VariantType.INS: QColor("#1e3a2e"),
    VariantType.INDEL: QColor("#3d3a1a"),
}


class VariantsPanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self._variants: list[Variant] = []
        self._filtered: list[Variant] = []
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        # Toolbar
        toolbar = QHBoxLayout()
        btn_open = QPushButton("Open VCF…")
        btn_open.setObjectName("primary")
        btn_open.clicked.connect(self._open_vcf)
        toolbar.addWidget(btn_open)

        toolbar.addWidget(QLabel("Filter:"))
        self._filter_input = QLineEdit()
        self._filter_input.setPlaceholderText("Gene, chromosome, or allele…")
        self._filter_input.textChanged.connect(self._apply_filter)
        toolbar.addWidget(self._filter_input, 1)

        self._status = QLabel("No VCF loaded.")
        self._status.setObjectName("subheading")
        toolbar.addWidget(self._status)
        layout.addLayout(toolbar)

        # Table
        self._table = QTableWidget()
        self._table.setColumnCount(8)
        self._table.setHorizontalHeaderLabels(
            ["CHROM", "POS", "REF", "ALT", "Type", "Gene", "Clinical Sig.", "gnomAD AF"]
        )
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self._table.horizontalHeader().setSectionResizeMode(6, QHeaderView.ResizeMode.Stretch)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self._table.currentCellChanged.connect(lambda row, *_: self._on_row_changed(row))
        layout.addWidget(self._table)

        # Detail panel
        detail_grp = QGroupBox("Variant Detail")
        detail_layout = QVBoxLayout(detail_grp)
        self._detail = QTextEdit()
        self._detail.setReadOnly(True)
        self._detail.setMaximumHeight(120)
        detail_layout.addWidget(self._detail)
        layout.addWidget(detail_grp)

    def _open_vcf(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Open VCF file", "", "VCF (*.vcf *.bcf);;All files (*)"
        )
        if not path:
            return
        try:
            with open(path) as f:
                self._variants = parse_vcf(f)
            self._filtered = self._variants[:]
            self._status.setText(f"{len(self._variants)} variants loaded.")
            self._populate_table(self._filtered)
        except Exception as e:
            QMessageBox.critical(self, "Parse error", str(e))

    def _apply_filter(self, text: str) -> None:
        t = text.lower()
        if not t:
            self._filtered = self._variants[:]
        else:
            self._filtered = [
                v for v in self._variants
                if t in v.chrom.lower()
                or t in v.gene.lower()
                or t in v.ref.lower()
                or any(t in a.lower() for a in v.alt)
                or (v.clinical_significance and t in v.clinical_significance.lower())
            ]
        self._status.setText(f"{len(self._filtered)}/{len(self._variants)} variants.")
        self._populate_table(self._filtered)

    def _populate_table(self, variants: list[Variant]) -> None:
        self._table.setRowCount(len(variants))
        SIG_COLORS = {
            "Pathogenic": QColor("#5c1111"),
            "Likely pathogenic": QColor("#5c3011"),
            "Likely benign": QColor("#1a2d3a"),
            "Benign": QColor("#1e3a2e"),
        }
        for row, v in enumerate(variants):
            row_color = TYPE_COLORS.get(v.variant_type, QColor("#1e293b"))
            vals = [
                v.chrom, str(v.pos), v.ref, ",".join(v.alt),
                v.variant_type.value, v.gene or "—",
                v.clinical_significance or "—",
                f"{v.gnomad_af:.2e}" if v.gnomad_af is not None else "—",
            ]
            for col, val in enumerate(vals):
                item = QTableWidgetItem(val)
                color = SIG_COLORS.get(v.clinical_significance or "", row_color) if col == 6 else row_color
                item.setBackground(color)
                self._table.setItem(row, col, item)

    def _on_row_changed(self, row: int) -> None:
        if row < 0 or row >= len(self._filtered):
            return
        v = self._filtered[row]
        lines = [
            f"ID: {v.id}",
            f"Position: {v.chrom}:{v.pos}",
            f"Ref → Alt: {v.ref} → {','.join(v.alt)}",
            f"Type: {v.variant_type.value}",
            f"Gene: {v.gene or '—'}",
            f"QUAL: {v.qual if v.qual is not None else '—'}",
            f"Filter: {';'.join(v.filter) or 'PASS'}",
            f"Clinical significance: {v.clinical_significance or '—'}",
            f"gnomAD AF: {v.gnomad_af if v.gnomad_af is not None else '—'}",
        ]
        if v.info:
            for k, val in list(v.info.items())[:8]:
                lines.append(f"{k}: {val}")
        self._detail.setPlainText("\n".join(lines))
