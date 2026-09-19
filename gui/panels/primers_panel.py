"""PCR primer design and restriction enzyme analysis panel."""

from __future__ import annotations
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import json
from pathlib import Path

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QTextEdit, QPushButton,
    QLabel, QSpinBox, QDoubleSpinBox, QGroupBox, QTableWidget,
    QTableWidgetItem, QProgressBar, QTabWidget, QMessageBox, QHeaderView,
    QComboBox, QInputDialog, QCheckBox, QSizePolicy,
)
from PyQt6.QtGui import QColor
from fastapi import HTTPException

from backend.routers.primers import (
    design_primers, restriction_sites, PrimerDesignRequest, RestrictionRequest,
    describe_tm_conditions, describe_pair_tm,
)
from backend.sequence_text import normalize_template
from ..table_header import fit_header_to_labels
from ..workers import Worker, worker_busy

# Presets stored in user config dir
_PRESETS_FILE = Path.home() / ".oligolia" / "primer_presets.json"
#: The combo's first row — a prompt, not a preset. Nothing is selected while
#: this is showing, which is what Delete has to know (#85.2).
_NO_PRESET = "— select preset —"
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
        self._last_digest = None  # most recent DigestResult, for restriction-ligation
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        # Template
        tmpl_grp = QGroupBox("Template Sequence")
        # Same defect as the CRISPR target box: ~370 px reserved for a field
        # capped at 70, with the primer table cramped underneath (#78.4b).
        tmpl_grp.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        tmpl_layout = QVBoxLayout(tmpl_grp)
        self._template = QTextEdit()
        self._template.setPlaceholderText("Paste or load template DNA sequence…")
        self._template.setMaximumHeight(70)
        tmpl_layout.addWidget(self._template)
        # Topology toggle — circular molecules need origin-spanning site search
        # (restriction sites / digest). Set automatically from a loaded
        # sequence, and user-adjustable for a pasted template.
        self._circular_check = QCheckBox("Circular topology")
        self._circular_check.setToolTip(
            "Search across the origin junction for restriction sites and digest "
            "(for plasmids and other circular molecules).")
        tmpl_layout.addWidget(self._circular_check)
        layout.addWidget(tmpl_grp)

        tabs = QTabWidget()

        # ── PCR Primer Design ────────────────────────────────────────────────
        pcr_widget = QWidget()
        pcr_layout = QVBoxLayout(pcr_widget)

        params_grp = QGroupBox("Parameters")
        # Six spin boxes and their six labels used to share one QHBoxLayout, so
        # at 1280 px the labels were the first thing to lose width and
        # "Product min (bp):" rendered as "Product min (bp" (#85.1). A grid wraps
        # them onto two rows of three pairs and gives each label column a minimum
        # taken from the label's own font metrics, so a label is never narrower
        # than the text it draws — on this platform's UI font or any other.
        params_layout = QGridLayout(params_grp)
        params_layout.setHorizontalSpacing(12)

        #: Labelled controls per row before wrapping to the next one.
        per_row = 3

        def param(index: int, label: str, widget: QWidget) -> QWidget:
            text = QLabel(label)
            text.setSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)
            row, col = divmod(index, per_row)
            params_layout.addWidget(text, row, col * 2)
            params_layout.addWidget(widget, row, col * 2 + 1)
            return widget

        def spin(lo: int, hi: int, val: int) -> QSpinBox:
            s = QSpinBox(); s.setRange(lo, hi); s.setValue(val)
            return s

        def dspin(lo: float, hi: float, val: float) -> QDoubleSpinBox:
            s = QDoubleSpinBox(); s.setRange(lo, hi); s.setValue(val)
            return s

        self._prod_min = param(0, "Product min (bp):", spin(50, 5000, 100))
        self._prod_max = param(1, "Product max (bp):", spin(50, 10000, 600))
        self._len_min = param(2, "Primer min:", spin(15, 35, 18))
        self._len_max = param(3, "Primer max:", spin(15, 35, 22))
        self._tm_min = param(4, "Tm min:", dspin(30, 90, 55.0))
        self._tm_max = param(5, "Tm max:", dspin(30, 90, 65.0))
        # The trailing stretch the old addStretch() gave the row: spare width
        # goes here rather than into the spin boxes.
        params_layout.setColumnStretch(per_row * 2, 1)
        pcr_layout.addWidget(params_grp)

        # Preset row
        preset_row = QHBoxLayout()
        preset_row.addWidget(QLabel("Preset:"))
        self._preset_combo = QComboBox()
        self._preset_combo.setMinimumWidth(180)
        self._preset_combo.currentTextChanged.connect(self._load_preset)
        self._preset_combo.currentTextChanged.connect(self._sync_preset_buttons)
        preset_row.addWidget(self._preset_combo)
        btn_save_preset = QPushButton("Save…")
        btn_save_preset.setObjectName("secondary")
        btn_save_preset.clicked.connect(self._save_preset)
        preset_row.addWidget(btn_save_preset)
        # Destructive, and until #85.2 offered against "— select preset —" and
        # against built-ins, neither of which it has ever been able to delete.
        # It keeps its styling for when it does have a target.
        self._btn_del_preset = QPushButton("Delete")
        self._btn_del_preset.setObjectName("secondary")
        self._btn_del_preset.clicked.connect(self._delete_preset)
        preset_row.addWidget(self._btn_del_preset)
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
        # "PRODUCT (BP)" was cut to "RODUCT (BP" at a flat 100 px per column —
        # the same header-truncation defect as the guide table, and it gets the
        # same platform-independent header (#78.3).
        fit_header_to_labels(self._pcr_table, stretch=(1, 2))
        self._pcr_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        pcr_layout.addWidget(self._pcr_table, 1)

        # A Tm is not a property of a sequence alone — the same primer reads 53.7
        # degC in this buffer and 64.0 degC in a PCR-like one. The table used to
        # show the number and say nothing about what it assumes (#79.4, #81).
        # Wording comes from the backend so it cannot drift from the number above it.
        self._pcr_tm_note = QLabel(describe_tm_conditions())
        self._pcr_tm_note.setObjectName("subheading")
        self._pcr_tm_note.setWordWrap(True)
        pcr_layout.addWidget(self._pcr_tm_note)
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

        # ── Assembly ──────────────────────────────────────────────────────────
        asm_widget = QWidget()
        asm_layout = QVBoxLayout(asm_widget)

        method_row = QHBoxLayout()
        method_row.addWidget(QLabel("Method:"))
        self._asm_method = QComboBox()
        self._asm_method.addItems(
            ["Gibson", "Golden Gate", "Restriction-ligation (from last digest)"])
        self._asm_method.currentTextChanged.connect(self._on_asm_method_changed)
        method_row.addWidget(self._asm_method, 1)

        self._asm_gg_enzyme = QComboBox()
        self._asm_gg_enzyme.addItems(["BsaI", "BsmBI", "Esp3I", "BbsI"])
        self._asm_gg_label = QLabel("Enzyme:")
        method_row.addWidget(self._asm_gg_label)
        method_row.addWidget(self._asm_gg_enzyme)

        self._asm_overlap_label = QLabel("Min overlap:")
        method_row.addWidget(self._asm_overlap_label)
        self._asm_overlap = QSpinBox(); self._asm_overlap.setRange(5, 60); self._asm_overlap.setValue(15)
        method_row.addWidget(self._asm_overlap)

        btn_asm = QPushButton("Assemble")
        btn_asm.setObjectName("primary")
        btn_asm.clicked.connect(self._run_assembly)
        method_row.addWidget(btn_asm)
        asm_layout.addLayout(method_row)

        self._asm_input = QTextEdit()
        self._asm_input.setPlaceholderText(
            "One fragment (Gibson) or part (Golden Gate) sequence per line…")
        self._asm_input.setMaximumHeight(110)
        asm_layout.addWidget(self._asm_input)

        self._asm_status = QLabel(""); self._asm_status.setObjectName("subheading")
        asm_layout.addWidget(self._asm_status)

        self._asm_result = QTextEdit()
        self._asm_result.setReadOnly(True)
        self._asm_result.setPlaceholderText("Assembled product appears here…")
        from PyQt6.QtGui import QFont
        self._asm_result.setFont(QFont("JetBrains Mono", 11))
        asm_layout.addWidget(self._asm_result)

        tabs.addTab(asm_widget, "Assembly")
        self._on_asm_method_changed(self._asm_method.currentText())

        layout.addWidget(tabs)

    def set_template(self, seq: str, is_circular: bool = False) -> None:
        self._template.setPlainText(seq)
        self._circular_check.setChecked(is_circular)

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
        self._preset_combo.addItem(_NO_PRESET)
        for name in self._all_presets():
            self._preset_combo.addItem(name)
        idx = self._preset_combo.findText(current)
        self._preset_combo.setCurrentIndex(max(0, idx))
        self._preset_combo.blockSignals(False)
        # Signals were blocked, so Delete has to be told by hand what the combo
        # now holds — including after a delete emptied the user's own presets.
        self._sync_preset_buttons()

    def _deletable_preset(self) -> str | None:
        """The selected preset if Delete could actually delete it, else None.

        The placeholder row is not a preset, and a built-in has always been
        refused with a warning box — so in both cases Delete has no target.
        """
        name = self._preset_combo.currentText()
        if not name or name == _NO_PRESET or name in _BUILTIN_PRESETS:
            return None
        return name if name in self._all_presets() else None

    def _sync_preset_buttons(self, *_: object) -> None:
        name = self._deletable_preset()
        self._btn_del_preset.setEnabled(name is not None)
        self._btn_del_preset.setToolTip(
            f"Delete the saved preset '{name}'." if name
            else "Select one of your own saved presets to delete it "
                 "(built-in presets cannot be deleted).")

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
        name = self._deletable_preset()
        if name is None:
            # The button is disabled in this state since #85.2; this stays as the
            # guard for any other caller.
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

    def _template_seq(self) -> str:
        """The template box read as one molecule.

        The same normalisation the routers apply, so the PCR, Digest and
        Restriction tabs all analyse the same sequence as each other and as the
        numbers they display — a GenBank ORIGIN paste or a CRLF sequence used to
        put the Restriction tab's positions and the Digest tab's cuts in two
        different coordinate systems (#86).
        """
        return normalize_template(self._template.toPlainText())

    def _run_pcr(self) -> None:
        template = self._template_seq()
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
        # Refuse rather than rebind a live QThread, which Qt aborts on (#54).
        if worker_busy(self, "_worker"):
            self._pcr_progress.hide()
            self._pcr_status.setText("Already designing primers — wait for that run to finish.")
            return
        self._worker = Worker(design_primers, req)
        self._worker.result.connect(self._on_pcr_done)
        self._worker.error.connect(lambda e: (self._pcr_progress.hide(), self._pcr_status.setText(f"Error: {e}")))
        self._worker.start()

    def _on_pcr_done(self, pairs: list) -> None:
        self._pcr_progress.hide()
        self._pcr_status.setText(f"{len(pairs)} primer pair{'s' if len(pairs) != 1 else ''} found.")
        # Say what these rows' Tm column assumes, taken from the pairs themselves
        # rather than from the defaults, so the caption stays true if the request
        # ever asks for a different buffer (#81).
        if pairs:
            self._pcr_tm_note.setText(describe_pair_tm(pairs[0]))
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
                    item.setFont(QFont("JetBrains Mono", 11))
                self._pcr_table.setItem(i, col, item)

    def _run_digest(self) -> None:
        template = self._template_seq()
        if not template:
            QMessageBox.warning(self, "No template", "Paste a template sequence first.")
            return
        enzyme_text = self._dig_enzymes.currentText().strip()
        if not enzyme_text:
            QMessageBox.warning(self, "No enzymes", "Enter one or more enzyme names.")
            return
        enzymes = [e.strip() for e in enzyme_text.replace("+", ",").split(",") if e.strip()]
        from backend.routers.primers import digest, DigestRequest, OVERHANG_AMBIGUOUS
        try:
            result = digest(DigestRequest(
                template=template, enzymes=enzymes,
                is_circular=self._circular_check.isChecked()))
            self._last_digest = result  # available to the Assembly tab
            n = len(result.fragments)
            cuts = len(result.cut_positions)
            # Two enzymes can cut the same bond and leave different ends (AvaI
            # and KpnI both cut pUC19's MCS at 412). The digest reports that end
            # as ambiguous instead of picking one (#88); say so here, because
            # the Assembly tab cannot ligate such an end and the table's
            # columns do not show end chemistry.
            shared = sorted({p for f in result.fragments
                             for p, t in ((f.start, f.left_overhang_type),
                                          (f.end, f.right_overhang_type))
                             if t == OVERHANG_AMBIGUOUS})
            note = ("  ⚠ ambiguous end at "
                    f"{', '.join(str(p) for p in shared)}: enzymes cutting the same "
                    "bond leave different ends" if shared else "")
            self._dig_status.setText(
                f"{cuts} cut site{'s' if cuts != 1 else ''} → {n} fragment{'s' if n != 1 else ''}  "
                f"(template: {result.template_length:,} bp)" + note
            )
            self._dig_table.setRowCount(n)
            for i, frag in enumerate(result.fragments):
                vals = [str(i + 1), f"{frag.length:,}", str(frag.start), str(frag.end)]
                for col, val in enumerate(vals):
                    item = QTableWidgetItem(val)
                    if i % 2 == 0:
                        item.setBackground(QColor("#1a2030"))
                    self._dig_table.setItem(i, col, item)
        except HTTPException as e:
            # The router's own refusal — an unknown enzyme, or a paste that is
            # not a nucleotide sequence (#87). Show its message, not
            # HTTPException's "400: …" stringification.
            self._dig_table.setRowCount(0)
            self._dig_status.setText(f"Error: {e.detail}")
        except Exception as e:
            self._dig_status.setText(f"Error: {e}")

    def _run_restriction(self) -> None:
        template = self._template_seq()
        if not template:
            return
        req = RestrictionRequest(template=template, is_circular=self._circular_check.isChecked())
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
        except HTTPException as e:
            # Same refusal the Digest tab gets for the same box, said the same
            # way, and with the table left empty rather than filled with
            # offsets into a FASTA header (#87).
            self._re_table.setRowCount(0)
            self._re_status.setText(f"Error: {e.detail}")
        except Exception as e:
            self._re_status.setText(f"Error: {e}")

    # ── Assembly ────────────────────────────────────────────────────────────

    def _on_asm_method_changed(self, method: str) -> None:
        is_gg = method.startswith("Golden")
        is_gibson = method.startswith("Gibson")
        is_ligation = method.startswith("Restriction")
        self._asm_gg_enzyme.setVisible(is_gg)
        self._asm_gg_label.setVisible(is_gg)
        self._asm_overlap.setVisible(is_gibson)
        self._asm_overlap_label.setVisible(is_gibson)
        self._asm_input.setVisible(not is_ligation)
        if is_ligation:
            self._asm_status.setText(
                "Ligates the fragments from the most recent Digest (in order) into a circle.")
        else:
            self._asm_status.setText("")

    def _asm_lines(self) -> list[str]:
        """One fragment per line, each normalised like any other template (#86).

        The newline is the separator here, so the lines are split first and the
        normalisation applied within each; a line that was only layout leaves
        nothing behind and is dropped, same as a blank one always was.
        """
        lines = (normalize_template(ln)
                 for ln in self._asm_input.toPlainText().splitlines())
        return [ln for ln in lines if ln]

    def _run_assembly(self) -> None:
        from fastapi import HTTPException
        method = self._asm_method.currentText()
        self._asm_result.clear()
        try:
            if method.startswith("Gibson"):
                from backend.routers.cloning import gibson, GibsonRequest
                frags = self._asm_lines()
                if len(frags) < 2:
                    QMessageBox.warning(self, "Need fragments", "Enter at least 2 fragment sequences (one per line).")
                    return
                res = gibson(GibsonRequest(fragments=frags, min_overlap=self._asm_overlap.value()))
            elif method.startswith("Golden"):
                from backend.routers.cloning import golden_gate, GoldenGateRequest
                parts = self._asm_lines()
                if len(parts) < 2:
                    QMessageBox.warning(self, "Need parts", "Enter at least 2 part sequences (one per line).")
                    return
                res = golden_gate(GoldenGateRequest(parts=parts, enzyme=self._asm_gg_enzyme.currentText()))
            else:  # Restriction-ligation from last digest
                from backend.routers.cloning import ligate, LigationRequest, LigationFragment
                if not self._last_digest or len(self._last_digest.fragments) < 2:
                    QMessageBox.warning(self, "No fragments",
                                        "Run a Digest producing ≥2 fragments first (Digest tab).")
                    return
                ordered = sorted(self._last_digest.fragments, key=lambda f: f.start)
                frags = [LigationFragment(**f.model_dump()) for f in ordered]
                # A linear digest has free ends on its terminal fragments and
                # cannot circularize; only close the loop when every end is a cut.
                free_ends = (ordered[0].left_overhang_type == "none"
                             or ordered[-1].right_overhang_type == "none")
                res = ligate(LigationRequest(fragments=frags, circular=not free_ends))
        except HTTPException as e:
            self._asm_status.setText(f"Assembly failed: {e.detail}")
            QMessageBox.critical(self, "Assembly failed", str(e.detail))
            return
        except Exception as e:
            self._asm_status.setText(f"Error: {e}")
            QMessageBox.critical(self, "Error", str(e))
            return

        product = res.product
        topo = "circular" if product.is_circular else "linear"
        njunc = len(res.junctions)
        self._asm_status.setText(
            f"✓ {topo} product · {product.length:,} bp · {njunc} junction{'s' if njunc != 1 else ''}"
            + (f"  ⚠ {len(res.warnings)} warning(s)" if res.warnings else "")
        )
        preview = product.seq if len(product.seq) <= 400 else product.seq[:400] + "…"
        lines = [f"// {product.name} ({topo}, {product.length:,} bp)", preview]
        if res.warnings:
            lines.append("")
            lines += [f"⚠ {w}" for w in res.warnings]
        self._asm_result.setPlainText("\n".join(lines))
