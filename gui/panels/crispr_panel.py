"""CRISPR guide RNA design panel."""

from __future__ import annotations
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, QPushButton,
    QTableWidget, QTableWidgetItem, QLabel, QComboBox, QSpinBox,
    QGroupBox, QProgressBar, QFileDialog, QMessageBox,
    QCheckBox, QSizePolicy,
)
from PyQt6.QtGui import QColor

from backend.models.crispr import CRISPRDesignRequest, CasType
from backend.routers.crispr import design_guides
from ..table_header import fit_header_to_labels, fitted_table
from ..workers import Worker, worker_busy


CAS_INFO = {
    "SpCas9": "NGG PAM · 20 nt guide · most widely used",
    "SpCas9-HF1": "NGG PAM · high fidelity · fewer off-targets",
    "AsCas12a": "TTTV PAM · 23 nt guide · staggered cuts",
    "LwaCas13a": "No PAM · RNA targeting · 22 nt guide",
}


class CRISPRPanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self._guides: list = []
        self._worker: Worker | None = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        # Target input
        target_grp = QGroupBox("Target Sequence")
        # The box held ~350 px for a field capped at 80: QTextEdit expands
        # vertically, so the group box claimed the panel's spare height and then
        # padded around a field that could not use it, while the results table
        # below stayed cramped. Fixed height hands that space to the table (#78.4).
        target_grp.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        target_layout = QVBoxLayout(target_grp)
        self._target_input = QTextEdit()
        self._target_input.setPlaceholderText("Paste target DNA sequence…")
        self._target_input.setMaximumHeight(80)
        target_layout.addWidget(self._target_input)
        target_layout.addWidget(QLabel("Tip: loaded sequence from the editor will be pre-filled."))
        layout.addWidget(target_grp)

        # Options
        opt_grp = QGroupBox("Options")
        opt_layout = QHBoxLayout(opt_grp)

        opt_layout.addWidget(QLabel("Cas type:"))
        self._cas_combo = QComboBox()
        for cas_val in ["SpCas9", "SpCas9-HF1", "AsCas12a", "LwaCas13a"]:
            self._cas_combo.addItem(cas_val)
        self._cas_combo.currentTextChanged.connect(
            lambda t: self._cas_info.setText(CAS_INFO.get(t, ""))
        )
        opt_layout.addWidget(self._cas_combo)

        opt_layout.addWidget(QLabel("Max guides:"))
        self._max_guides = QSpinBox(); self._max_guides.setRange(1, 50); self._max_guides.setValue(10)
        opt_layout.addWidget(self._max_guides)

        self._check_off = QCheckBox("Check off-targets")
        self._check_off.setToolTip(
            "Scan the target sequence for near-match sites and score guide "
            "specificity (MIT). Not available for LwaCas13a (RNA target)."
        )
        opt_layout.addWidget(self._check_off)

        opt_layout.addStretch()
        layout.addWidget(opt_grp)

        self._cas_info = QLabel(CAS_INFO["SpCas9"])
        self._cas_info.setObjectName("subheading")
        layout.addWidget(self._cas_info)

        # Run button + progress
        btn_row = QHBoxLayout()
        self._btn_design = QPushButton("Design Guides")
        self._btn_design.setObjectName("primary")
        self._btn_design.clicked.connect(self._run_design)
        btn_row.addWidget(self._btn_design)
        self._btn_export = QPushButton("Export TSV")
        self._btn_export.clicked.connect(self._export_tsv)
        self._btn_export.setEnabled(False)
        btn_row.addWidget(self._btn_export)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        self._progress = QProgressBar()
        self._progress.setRange(0, 0)  # indeterminate
        self._progress.hide()
        layout.addWidget(self._progress)

        self._status = QLabel("")
        self._status.setObjectName("subheading")
        layout.addWidget(self._status)

        # Results table
        self._table = fitted_table()
        self._table.setColumnCount(9)
        self._table.setHorizontalHeaderLabels(
            ["#", "Guide Sequence (5'→3')", "PAM", "Position", "Strand", "GC%",
             "On-target", "Off-targets", "Specificity"]
        )
        # Every column was a flat 100 px, so "ON-TARGET SCORE" rendered as
        # "N-TARGET SCOR" — one of the two numbers a guide is picked on (#78.3).
        # Sizing to contents measures the header text too, but in Qt's choice of
        # font and case, which fit here and clipped by 3 px on macOS; the header
        # below raises every section to what its own upper-cased label needs in
        # the live fontMetrics, so it holds on whatever platform runs it.
        fit_header_to_labels(self._table, stretch=(1,))
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        layout.addWidget(self._table, 1)

    #: Above this, say so instead of quietly shortening the search space.
    LONG_TARGET_NT = 2000

    def set_target(self, seq: str) -> None:
        """Pre-fill target from the loaded sequence.

        This used to do `seq[:2000]` "for display", which silently dropped every
        site beyond 2 kb from the search and told the user nothing — a 5 kb insert
        was designed against its first 2 kb (#64.4). The whole sequence is loaded
        now; long ones get a note on the status line rather than a haircut.
        """
        self._target_input.setPlainText(seq)
        if len(seq) > self.LONG_TARGET_NT:
            self._status.setText(
                f"Target is {len(seq):,} nt — the whole sequence will be searched. "
                "Trim the box if you only want part of it."
            )
        else:
            self._status.setText("")

    def _run_design(self) -> None:
        target = self._target_input.toPlainText().strip().upper().replace(" ", "").replace("\n", "")
        if not target:
            QMessageBox.warning(self, "No target", "Paste a target DNA sequence first.")
            return
        self._btn_design.setEnabled(False)
        self._progress.show()
        self._table.setRowCount(0)
        self._guides.clear()
        self._btn_export.setEnabled(False)

        cas_str = self._cas_combo.currentText()
        # guide_length is deliberately NOT passed. The panel has no control for
        # it, and the router uses each nuclease's canonical length (20 for Cas9,
        # 23 for Cas12a, 22 for Cas13) unless a caller explicitly asks for one.
        # Passing 20 here would have silently shortened every Cas12a guide (#64).
        req = CRISPRDesignRequest(
            target_sequence=target,
            cas_type=CasType(cas_str),
            max_guides=self._max_guides.value(),
            check_off_targets=self._check_off.isChecked(),
        )
        # Refuse rather than rebind a live QThread, which Qt aborts on (#54).
        if worker_busy(self, "_worker"):
            self._progress.hide()
            self._btn_design.setEnabled(True)
            self._status.setText("Already designing guides — wait for that run to finish.")
            return
        self._worker = Worker(design_guides, req)
        self._worker.result.connect(self._on_design_done)
        self._worker.error.connect(self._on_design_error)
        self._worker.start()

    def _on_design_done(self, response) -> None:
        self._progress.hide()
        self._btn_design.setEnabled(True)
        self._guides = response.guides
        self._status.setText(
            f"{len(response.guides)} guides shown (of {response.total_candidates} candidates) "
            f"for {response.target_length} bp target"
        )
        self._table.setRowCount(len(response.guides))

        for i, g in enumerate(response.guides):
            score = g.on_target_score or 0.0
            color = QColor("#1e3a2e") if score >= 0.6 else QColor("#3d2a00") if score >= 0.4 else QColor("#3d1a1a")

            if g.off_target_summary is not None:
                s = g.off_target_summary
                off_text = (
                    f"{s.get('0', 0)} exact · {s.get('1', 0)} (1mm) · {s.get('2', 0)} (2mm)"
                )
                spec_text = f"{g.specificity_score:.0f}" if g.specificity_score is not None else "—"
            else:
                off_text = "—"
                spec_text = "—"

            for col, val in enumerate([
                str(i + 1), g.sequence, g.pam,
                str(g.position), g.strand,
                f"{g.gc_content:.0f}%",
                f"{score:.2f}",
                off_text,
                spec_text,
            ]):
                item = QTableWidgetItem(val)
                item.setBackground(color)
                if col == 1:
                    from PyQt6.QtGui import QFont
                    item.setFont(QFont("JetBrains Mono", 11))
                self._table.setItem(i, col, item)

        self._btn_export.setEnabled(True)

    def _on_design_error(self, err: str) -> None:
        self._progress.hide()
        self._btn_design.setEnabled(True)
        self._status.setText(f"Error: {err}")

    def _export_tsv(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Export guides", "", "TSV (*.tsv);;All files (*)")
        if not path:
            return
        with open(path, "w") as f:
            f.write(
                "sequence\tpam\tposition\tstrand\tgc_content\ton_target_score\t"
                "off_target_exact\toff_target_1mm\toff_target_2mm\toff_target_3mm\tspecificity\n"
            )
            for g in self._guides:
                s = g.off_target_summary or {}
                f.write(
                    f"{g.sequence}\t{g.pam}\t{g.position}\t{g.strand}\t{g.gc_content}\t"
                    f"{g.on_target_score or ''}\t"
                    f"{s.get('0', '')}\t{s.get('1', '')}\t{s.get('2', '')}\t{s.get('3', '')}\t"
                    f"{g.specificity_score if g.specificity_score is not None else ''}\n"
                )
