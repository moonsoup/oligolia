"""Workflow builder panel — chain operations via the backend engine (issue #10).

This wires the *existing* backend/workflow engine into the GUI: build an
ordered pipeline of steps, seed it with a sequence, run it, and see per-step
status/results, plus save/load as .ogo.

Scope note: the "order" step and any vendor/live-ordering path are
deliberately NOT exposed here — that touches real money / external accounts
and needs explicit direction first. Only the offline computational steps the
engine supports are offered.
"""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QComboBox,
    QTextEdit, QGroupBox, QTableWidget, QTableWidgetItem, QHeaderView,
    QMessageBox, QFileDialog, QLineEdit,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor

from backend.workflow import (
    StepStatus, StepType, Workflow, WorkflowStep, run_workflow, save_ogo, load_ogo,
)

# Steps exposed in the GUI. "order"/db_search/codon_optimize/msa are excluded:
# order/vendor pathways are out of scope for now, and the others aren't wired
# in the engine yet (they report a clear error).
_GUI_STEPS = [
    StepType.CRISPR_DESIGN,
    StepType.OFF_TARGET,
    StepType.PRIMER_DESIGN,
    StepType.RESTRICTION_DIGEST,
    StepType.EXPORT,
]

# Sensible starting params per step so users aren't writing JSON from scratch.
_DEFAULT_PARAMS = {
    StepType.CRISPR_DESIGN: {"max_guides": 5, "check_off_targets": False},
    StepType.OFF_TARGET: {"min_specificity": 0},
    StepType.PRIMER_DESIGN: {},
    StepType.RESTRICTION_DIGEST: {"enzymes": ["EcoRI", "BamHI"]},
    StepType.EXPORT: {"source": "guides"},
}

_STATUS_COLOR = {
    StepStatus.COMPLETE: "#1e3a2e",
    StepStatus.FAILED: "#3d1a1a",
    StepStatus.SKIPPED: "#2a2a2a",
    StepStatus.RUNNING: "#3d2a00",
    StepStatus.PENDING: "#1a2030",
}


def _result_summary(step: WorkflowStep) -> str:
    """Short human summary of a finished step's result/error."""
    if step.status == StepStatus.FAILED:
        return step.error or "failed"
    if step.status != StepStatus.COMPLETE:
        return ""
    r = step.result
    if not isinstance(r, dict):
        return "done"
    if step.type == StepType.CRISPR_DESIGN:
        return f"{len(r.get('guides', []))} guides ({r.get('total_candidates', 0)} candidates)"
    if step.type == StepType.OFF_TARGET:
        return f"{r.get('passed', 0)}/{r.get('scanned', 0)} guides pass"
    if step.type == StepType.PRIMER_DESIGN:
        return f"{r.get('count', 0)} primer pairs"
    if step.type == StepType.RESTRICTION_DIGEST:
        return f"{len(r.get('fragments', []))} fragments"
    if step.type == StepType.EXPORT:
        return r.get("filename", "exported")
    return "done"


class WorkflowPanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self._last_wf: Workflow | None = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        # ── Seed sequence ─────────────────────────────────────────────────
        seq_grp = QGroupBox("Input Sequence (seeds the workflow context)")
        seq_layout = QVBoxLayout(seq_grp)
        name_row = QHBoxLayout()
        name_row.addWidget(QLabel("Workflow name:"))
        self._name_edit = QLineEdit("My workflow")
        name_row.addWidget(self._name_edit, 1)
        seq_layout.addLayout(name_row)
        self._seq_input = QTextEdit()
        self._seq_input.setPlaceholderText("Paste/load a DNA sequence, or select one in the Sequences tab…")
        self._seq_input.setMaximumHeight(60)
        seq_layout.addWidget(self._seq_input)
        layout.addWidget(seq_grp)

        # ── Steps builder ─────────────────────────────────────────────────
        steps_grp = QGroupBox("Steps")
        steps_layout = QVBoxLayout(steps_grp)

        add_row = QHBoxLayout()
        self._step_combo = QComboBox()
        for st in _GUI_STEPS:
            self._step_combo.addItem(st.value, st)
        add_row.addWidget(self._step_combo, 1)
        btn_add = QPushButton("Add step")
        btn_add.setObjectName("primary")
        btn_add.clicked.connect(self._add_step)
        add_row.addWidget(btn_add)
        btn_up = QPushButton("Move ↑"); btn_up.clicked.connect(lambda: self._move_step(-1))
        btn_down = QPushButton("Move ↓"); btn_down.clicked.connect(lambda: self._move_step(1))
        btn_rm = QPushButton("Remove"); btn_rm.setObjectName("danger"); btn_rm.clicked.connect(self._remove_step)
        add_row.addWidget(btn_up); add_row.addWidget(btn_down); add_row.addWidget(btn_rm)
        steps_layout.addLayout(add_row)

        self._steps_table = QTableWidget()
        self._steps_table.setColumnCount(2)
        self._steps_table.setHorizontalHeaderLabels(["Step", "Parameters (JSON)"])
        self._steps_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self._steps_table.setMaximumHeight(150)
        steps_layout.addWidget(self._steps_table)
        layout.addWidget(steps_grp)

        # ── Run / save / load ─────────────────────────────────────────────
        run_row = QHBoxLayout()
        btn_run = QPushButton("Run Workflow")
        btn_run.setObjectName("primary")
        btn_run.clicked.connect(self._run)
        run_row.addWidget(btn_run)
        btn_save = QPushButton("Save .ogo…"); btn_save.clicked.connect(self._save_ogo)
        btn_load = QPushButton("Load .ogo…"); btn_load.clicked.connect(self._load_ogo)
        run_row.addWidget(btn_save); run_row.addWidget(btn_load)
        run_row.addStretch()
        self._status = QLabel(""); self._status.setObjectName("subheading")
        run_row.addWidget(self._status)
        layout.addLayout(run_row)

        # ── Results ───────────────────────────────────────────────────────
        self._results_table = QTableWidget()
        self._results_table.setColumnCount(3)
        self._results_table.setHorizontalHeaderLabels(["Step", "Status", "Detail"])
        self._results_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self._results_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._results_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self._results_table.itemSelectionChanged.connect(self._show_result_detail)
        layout.addWidget(self._results_table)

        self._detail = QTextEdit()
        self._detail.setReadOnly(True)
        self._detail.setPlaceholderText("Select a completed step to inspect its full result…")
        self._detail.setMaximumHeight(140)
        layout.addWidget(self._detail)

    # ── Cross-panel hook ──────────────────────────────────────────────────
    def set_sequence(self, seq: str) -> None:
        self._seq_input.setPlainText(seq)

    # ── Steps table editing ───────────────────────────────────────────────
    def _add_step(self, _checked: bool = False, step_type: StepType | None = None,
                  params: dict | None = None) -> None:
        st = step_type or self._step_combo.currentData()
        row = self._steps_table.rowCount()
        self._steps_table.insertRow(row)
        type_item = QTableWidgetItem(st.value)
        type_item.setFlags(type_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
        self._steps_table.setItem(row, 0, type_item)
        params = _DEFAULT_PARAMS.get(st, {}) if params is None else params
        self._steps_table.setItem(row, 1, QTableWidgetItem(json.dumps(params)))

    def _move_step(self, delta: int) -> None:
        row = self._steps_table.currentRow()
        new = row + delta
        if row < 0 or not (0 <= new < self._steps_table.rowCount()):
            return
        for col in range(2):
            a = self._steps_table.takeItem(row, col)
            b = self._steps_table.takeItem(new, col)
            self._steps_table.setItem(row, col, b)
            self._steps_table.setItem(new, col, a)
        self._steps_table.setCurrentCell(new, 0)

    def _remove_step(self) -> None:
        row = self._steps_table.currentRow()
        if row >= 0:
            self._steps_table.removeRow(row)

    def _build_workflow(self) -> Workflow:
        steps: list[WorkflowStep] = []
        for row in range(self._steps_table.rowCount()):
            st_value = self._steps_table.item(row, 0).text()
            raw = self._steps_table.item(row, 1).text().strip() or "{}"
            params = json.loads(raw)  # may raise; caller handles
            steps.append(WorkflowStep(id=f"s{row + 1}", type=StepType(st_value), params=params))
        return Workflow(name=self._name_edit.text() or "workflow", steps=steps)

    # ── Run ────────────────────────────────────────────────────────────────
    def _run(self) -> None:
        if self._steps_table.rowCount() == 0:
            QMessageBox.warning(self, "No steps", "Add at least one step.")
            return
        try:
            wf = self._build_workflow()
        except json.JSONDecodeError as e:
            QMessageBox.critical(self, "Invalid parameters", f"A step's JSON parameters are invalid:\n{e}")
            return

        seq = self._seq_input.toPlainText().strip().upper().replace(" ", "").replace("\n", "")
        context = {"sequence": seq} if seq else {}
        run_workflow(wf, context)
        self._last_wf = wf
        self._render_results(wf)
        done = sum(1 for s in wf.steps if s.status == StepStatus.COMPLETE)
        self._status.setText(f"{done}/{len(wf.steps)} steps completed"
                             + ("  ·  workflow complete" if wf.is_complete else ""))

    def _render_results(self, wf: Workflow) -> None:
        self._results_table.setRowCount(len(wf.steps))
        for i, step in enumerate(wf.steps):
            name = QTableWidgetItem(f"{i + 1}. {step.type.value}")
            status = QTableWidgetItem(step.status.value)
            status.setBackground(QColor(_STATUS_COLOR.get(step.status, "#1a2030")))
            detail = QTableWidgetItem(_result_summary(step))
            for col, item in enumerate([name, status, detail]):
                self._results_table.setItem(i, col, item)

    def _show_result_detail(self) -> None:
        if not self._last_wf:
            return
        rows = self._results_table.selectionModel().selectedRows()
        if not rows:
            return
        step = self._last_wf.steps[rows[0].row()]
        if step.status == StepStatus.FAILED:
            self._detail.setPlainText(f"// {step.type.value} — FAILED\n{step.error}")
        else:
            self._detail.setPlainText(json.dumps(step.result, indent=2, default=str))

    # ── .ogo save/load ──────────────────────────────────────────────────────
    def _save_ogo(self) -> None:
        try:
            wf = self._last_wf or self._build_workflow()
        except json.JSONDecodeError as e:
            QMessageBox.critical(self, "Invalid parameters", str(e))
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save workflow", f"{wf.name}.ogo", "Oligolia workflow (*.ogo)")
        if not path:
            return
        save_ogo(wf, path)
        self._status.setText(f"Saved {os.path.basename(path)}")

    def _load_ogo(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Load workflow", "", "Oligolia workflow (*.ogo);;All files (*)")
        if not path:
            return
        try:
            wf = load_ogo(path)
        except Exception as e:
            QMessageBox.critical(self, "Load failed", str(e))
            return
        self._name_edit.setText(wf.name)
        self._steps_table.setRowCount(0)
        for step in wf.steps:
            try:
                st = StepType(step.type)
            except ValueError:
                continue
            self._add_step(step_type=st, params=step.params)
        self._last_wf = wf
        self._render_results(wf)
        self._status.setText(f"Loaded {os.path.basename(path)} ({len(wf.steps)} steps)")
