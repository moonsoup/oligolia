"""Sequence viewer and editor panel."""

from __future__ import annotations
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QListWidget, QListWidgetItem,
    QTextEdit, QLabel, QPushButton, QComboBox, QLineEdit, QGroupBox,
    QFormLayout, QSpinBox, QMessageBox, QFileDialog, QApplication,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QFont, QColor, QTextCharFormat, QSyntaxHighlighter

from Bio.Seq import Seq
from backend.models.sequence import Sequence, MoleculeType
from backend.formats import read_fasta, read_fastq, read_genbank
import re


class DNAHighlighter(QSyntaxHighlighter):
    """Color nucleotide bases: A=green, T=red, G=gold, C=blue."""
    COLORS = {"A": "#4ade80", "T": "#f87171", "G": "#fbbf24", "C": "#60a5fa",
               "U": "#c084fc", "-": "#475569"}

    def highlightBlock(self, text: str) -> None:
        for i, ch in enumerate(text.upper()):
            color = self.COLORS.get(ch)
            if color:
                fmt = QTextCharFormat()
                fmt.setForeground(QColor(color))
                self.setFormat(i, 1, fmt)


class SequencePanel(QWidget):
    sequence_selected = pyqtSignal(object)  # emits Sequence

    def __init__(self) -> None:
        super().__init__()
        self._sequences: dict[str, Sequence] = {}
        self._active: Sequence | None = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ── Left: sequence list ──────────────────────────────────────────────
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(8, 8, 8, 8)

        title = QLabel("Sequences")
        title.setObjectName("heading")
        left_layout.addWidget(title)

        self._list = QListWidget()
        self._list.currentItemChanged.connect(self._on_seq_selected)
        left_layout.addWidget(self._list)

        btn_row = QHBoxLayout()
        btn_open = QPushButton("Open file…")
        btn_open.clicked.connect(self._open_file)
        btn_del = QPushButton("Remove")
        btn_del.setObjectName("danger")
        btn_del.clicked.connect(self._remove_selected)
        btn_row.addWidget(btn_open)
        btn_row.addWidget(btn_del)
        left_layout.addLayout(btn_row)

        # Paste / new sequence
        paste_grp = QGroupBox("Add sequence manually")
        paste_layout = QFormLayout(paste_grp)
        self._paste_id = QLineEdit()
        self._paste_id.setPlaceholderText("e.g. MyGene")
        self._paste_seq = QTextEdit()
        self._paste_seq.setPlaceholderText("Paste sequence here…")
        self._paste_seq.setMaximumHeight(60)
        self._paste_mol = QComboBox()
        self._paste_mol.addItems(["DNA", "RNA", "PROTEIN"])
        btn_paste = QPushButton("Add")
        btn_paste.setObjectName("primary")
        btn_paste.clicked.connect(self._add_manual)
        paste_layout.addRow("ID:", self._paste_id)
        paste_layout.addRow("Type:", self._paste_mol)
        paste_layout.addRow("Seq:", self._paste_seq)
        paste_layout.addRow(btn_paste)
        left_layout.addWidget(paste_grp)

        splitter.addWidget(left)

        # ── Right: viewer + editor ───────────────────────────────────────────
        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(8, 8, 8, 8)

        # Sequence info bar
        self._info_label = QLabel("No sequence selected")
        self._info_label.setObjectName("subheading")
        right_layout.addWidget(self._info_label)

        # Sequence display
        self._seq_display = QTextEdit()
        self._seq_display.setReadOnly(True)
        self._seq_display.setFont(QFont("JetBrains Mono,Fira Code,Courier New", 11))
        self._highlighter = DNAHighlighter(self._seq_display.document())
        right_layout.addWidget(self._seq_display)

        # Edit operations
        edit_grp = QGroupBox("Edit Operation")
        edit_layout = QVBoxLayout(edit_grp)

        op_row = QHBoxLayout()
        self._op_combo = QComboBox()
        for op in ["reverse_complement", "complement", "translate", "transcribe",
                   "back_transcribe", "insert", "delete", "replace"]:
            self._op_combo.addItem(op.replace("_", " ").title(), op)
        self._op_combo.currentIndexChanged.connect(self._update_op_fields)
        op_row.addWidget(self._op_combo, 1)
        btn_apply = QPushButton("Apply")
        btn_apply.setObjectName("primary")
        btn_apply.clicked.connect(self._apply_op)
        op_row.addWidget(btn_apply)
        edit_layout.addLayout(op_row)

        # Position fields (shown/hidden by operation)
        self._pos_widget = QWidget()
        pos_layout = QHBoxLayout(self._pos_widget)
        pos_layout.setContentsMargins(0, 0, 0, 0)
        pos_layout.addWidget(QLabel("Start:"))
        self._pos_start = QSpinBox(); self._pos_start.setMaximum(999999)
        pos_layout.addWidget(self._pos_start)
        pos_layout.addWidget(QLabel("End:"))
        self._pos_end = QSpinBox(); self._pos_end.setMaximum(999999)
        pos_layout.addWidget(self._pos_end)
        self._pos_widget.hide()
        edit_layout.addWidget(self._pos_widget)

        self._insert_widget = QWidget()
        ins_layout = QHBoxLayout(self._insert_widget)
        ins_layout.setContentsMargins(0, 0, 0, 0)
        ins_layout.addWidget(QLabel("Sequence:"))
        self._insert_seq = QLineEdit()
        self._insert_seq.setPlaceholderText("Bases to insert/replace")
        ins_layout.addWidget(self._insert_seq)
        self._insert_widget.hide()
        edit_layout.addWidget(self._insert_widget)

        right_layout.addWidget(edit_grp)

        # Motif search
        motif_grp = QGroupBox("Find Motif (IUPAC)")
        motif_layout = QHBoxLayout(motif_grp)
        self._motif_input = QLineEdit()
        self._motif_input.setPlaceholderText("e.g. GAATTC or RRNYY")
        motif_layout.addWidget(self._motif_input)
        btn_motif = QPushButton("Search")
        btn_motif.clicked.connect(self._find_motif)
        motif_layout.addWidget(btn_motif)
        self._motif_result = QLabel("")
        motif_layout.addWidget(self._motif_result)
        right_layout.addWidget(motif_grp)

        # Result display
        self._result_display = QTextEdit()
        self._result_display.setReadOnly(True)
        self._result_display.setPlaceholderText("Results will appear here…")
        self._result_display.setMaximumHeight(120)
        right_layout.addWidget(self._result_display)

        result_btns = QHBoxLayout()
        btn_copy = QPushButton("Copy result")
        btn_copy.clicked.connect(lambda: QApplication.clipboard().setText(self._result_display.toPlainText()))
        btn_save_fasta = QPushButton("Save result as FASTA")
        btn_save_fasta.clicked.connect(self._save_result_fasta)
        result_btns.addWidget(btn_copy)
        result_btns.addWidget(btn_save_fasta)
        result_btns.addStretch()
        right_layout.addLayout(result_btns)

        splitter.addWidget(right)
        splitter.setSizes([280, 700])
        layout.addWidget(splitter)

    def _update_op_fields(self) -> None:
        op = self._op_combo.currentData()
        self._pos_widget.setVisible(op in ("delete", "replace"))
        self._insert_widget.setVisible(op in ("insert", "replace"))

    def add_sequence(self, seq: Sequence) -> None:
        self._sequences[seq.id] = seq
        item = QListWidgetItem(f"{'🔵' if seq.molecule_type == MoleculeType.DNA else '🟡' if seq.molecule_type == MoleculeType.RNA else '🟣'} {seq.id}")
        item.setData(Qt.ItemDataRole.UserRole, seq.id)
        item.setToolTip(f"{seq.length:,} bp · {seq.molecule_type.value}")
        self._list.addItem(item)
        self._list.setCurrentItem(item)

    def _on_seq_selected(self, item: QListWidgetItem | None, _: object = None) -> None:
        if not item:
            return
        seq_id = item.data(Qt.ItemDataRole.UserRole)
        seq = self._sequences.get(seq_id)
        if not seq:
            return
        self._active = seq
        self._info_label.setText(
            f"{seq.name or seq.id}  ·  {seq.molecule_type.value}  ·  {seq.length:,} bp  ·  "
            f"GC: {self._gc(seq.seq):.1f}%"
        )
        self._seq_display.setPlainText(seq.seq)
        self.sequence_selected.emit(seq)

    def _gc(self, seq: str) -> float:
        s = seq.upper()
        return (s.count("G") + s.count("C")) / len(s) * 100 if s else 0.0

    def _apply_op(self) -> None:
        if not self._active:
            QMessageBox.warning(self, "No sequence", "Select a sequence first.")
            return
        op = self._op_combo.currentData()
        seq_str = self._active.seq

        try:
            bio = Seq(seq_str)
            if op == "reverse_complement":
                result = str(bio.reverse_complement())
                msg = "Reverse complement"
            elif op == "complement":
                result = str(bio.complement())
                msg = "Complement"
            elif op == "translate":
                result = str(bio.translate(to_stop=True))
                msg = f"Translation ({len(seq_str)} nt → {len(result)} aa)"
            elif op == "transcribe":
                result = str(bio.transcribe())
                msg = "Transcription (DNA → RNA)"
            elif op == "back_transcribe":
                result = str(bio.back_transcribe())
                msg = "Back-transcription (RNA → DNA)"
            elif op == "insert":
                pos = self._pos_start.value()
                ins = self._insert_seq.text().upper()
                result = seq_str[:pos] + ins + seq_str[pos:]
                msg = f"Inserted {len(ins)} bases at position {pos}"
            elif op == "delete":
                start, end = self._pos_start.value(), self._pos_end.value()
                result = seq_str[:start] + seq_str[end:]
                msg = f"Deleted [{start}:{end}]"
            elif op == "replace":
                start, end = self._pos_start.value(), self._pos_end.value()
                rep = self._insert_seq.text().upper()
                result = seq_str[:start] + rep + seq_str[end:]
                msg = f"Replaced [{start}:{end}] → {len(rep)} bases"
            else:
                return

            self._result_display.setPlainText(f"// {msg}\n{result}")
            self._last_result = result
            self._last_result_id = f"{self._active.id}_{op}"

        except Exception as e:
            QMessageBox.critical(self, "Operation failed", str(e))

    def _find_motif(self) -> None:
        if not self._active:
            return
        motif = self._motif_input.text().upper()
        if not motif:
            return
        IUPAC = {"R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
                 "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
                 "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
        pattern = "".join(IUPAC.get(c, re.escape(c)) for c in motif)
        positions = [m.start() for m in re.finditer(f"(?={pattern})", self._active.seq.upper())]
        self._motif_result.setText(f"{len(positions)} hit{'s' if len(positions) != 1 else ''}")
        pos_str = ", ".join(str(p) for p in positions[:20])
        if positions:
            pos_str += ("…" if len(positions) > 20 else "")
            self._result_display.setPlainText(f"// Motif '{motif}' — {len(positions)} occurrences\nPositions: {pos_str}")

    def _save_result_fasta(self) -> None:
        result = self._result_display.toPlainText()
        if not result or result.startswith("//"):
            lines = result.splitlines()
            seq_lines = [ln for ln in lines if ln and not ln.startswith("//")]
            seq_text = "".join(seq_lines)
        else:
            seq_text = result
        if not seq_text:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save FASTA", "", "FASTA (*.fa *.fasta);;All files (*)")
        if path:
            seq_id = getattr(self, "_last_result_id", "result")
            with open(path, "w") as f:
                f.write(f">{seq_id}\n{seq_text}\n")

    def _open_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Open sequence file", "",
            "All bioinformatics (*.fasta *.fa *.fna *.faa *.fastq *.fq *.gb *.gbk *.embl);;"
            "FASTA (*.fasta *.fa *.fna *.faa);;"
            "FASTQ (*.fastq *.fq);;"
            "GenBank (*.gb *.gbk);;"
            "EMBL (*.embl);;All files (*)"
        )
        if not path:
            return
        ext = path.rsplit(".", 1)[-1].lower()
        try:
            readers = {"fasta": read_fasta, "fa": read_fasta, "fna": read_fasta,
                       "faa": read_fasta, "fastq": read_fastq, "fq": read_fastq,
                       "gb": read_genbank, "gbk": read_genbank, "embl": read_genbank}
            reader = readers.get(ext, read_fasta)
            with open(path) as f:
                seqs = reader(f)
            for seq in seqs:
                self.add_sequence(seq)
        except Exception as e:
            QMessageBox.critical(self, "Import failed", str(e))

    def _remove_selected(self) -> None:
        item = self._list.currentItem()
        if item:
            seq_id = item.data(Qt.ItemDataRole.UserRole)
            self._sequences.pop(seq_id, None)
            self._list.takeItem(self._list.row(item))

    def _add_manual(self) -> None:
        seq_id = self._paste_id.text().strip() or "manual"
        raw = self._paste_seq.toPlainText().strip().upper().replace(" ", "").replace("\n", "")
        mol = MoleculeType(self._paste_mol.currentText())
        if not raw:
            QMessageBox.warning(self, "Empty", "Paste a sequence first.")
            return
        seq = Sequence(id=seq_id, name=seq_id, seq=raw, molecule_type=mol, description="")
        self.add_sequence(seq)
        self._paste_id.clear()
        self._paste_seq.clear()

    @property
    def sequences(self) -> dict[str, Sequence]:
        return self._sequences
