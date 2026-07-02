"""Sequence alignment panel — pairwise and MSA."""

from __future__ import annotations
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, QPushButton,
    QLabel, QComboBox, QTabWidget, QTableWidget, QTableWidgetItem,
    QProgressBar, QMessageBox,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QColor
from Bio import Align
from ..workers import Worker


class AlignmentPanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self._msa_worker: Worker | None = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        tabs = QTabWidget()

        # ── Pairwise ──────────────────────────────────────────────────────────
        pair_widget = QWidget()
        pair_layout = QVBoxLayout(pair_widget)

        mode_row = QHBoxLayout()
        self._mode = QComboBox()
        self._mode.addItem("Global (Needleman-Wunsch)", "global")
        self._mode.addItem("Local (Smith-Waterman)", "local")
        mode_row.addWidget(QLabel("Mode:"))
        mode_row.addWidget(self._mode)
        mode_row.addStretch()
        pair_layout.addLayout(mode_row)

        self._seq1 = QTextEdit()
        self._seq1.setPlaceholderText("Sequence 1…")
        self._seq1.setMaximumHeight(70)
        self._seq1.setFont(QFont("JetBrains Mono", 11))
        pair_layout.addWidget(QLabel("Sequence 1:"))
        pair_layout.addWidget(self._seq1)

        self._seq2 = QTextEdit()
        self._seq2.setPlaceholderText("Sequence 2…")
        self._seq2.setMaximumHeight(70)
        self._seq2.setFont(QFont("JetBrains Mono", 11))
        pair_layout.addWidget(QLabel("Sequence 2:"))
        pair_layout.addWidget(self._seq2)

        btn_pair = QPushButton("Align")
        btn_pair.setObjectName("primary")
        btn_pair.clicked.connect(self._run_pairwise)
        pair_layout.addWidget(btn_pair)

        # Stats row
        stats_row = QHBoxLayout()
        self._stat_score = QLabel("Score: —")
        self._stat_identity = QLabel("Identity: —")
        self._stat_gaps = QLabel("Gaps: —")
        self._stat_len = QLabel("Length: —")
        for lbl in [self._stat_score, self._stat_identity, self._stat_gaps, self._stat_len]:
            lbl.setObjectName("subheading")
            stats_row.addWidget(lbl)
        stats_row.addStretch()
        pair_layout.addLayout(stats_row)

        self._pair_result = QTextEdit()
        self._pair_result.setReadOnly(True)
        self._pair_result.setFont(QFont("JetBrains Mono", 11))
        pair_layout.addWidget(self._pair_result)

        tabs.addTab(pair_widget, "Pairwise")

        # ── MSA ───────────────────────────────────────────────────────────────
        msa_widget = QWidget()
        msa_layout = QVBoxLayout(msa_widget)

        msa_layout.addWidget(QLabel("Enter sequences (one per line, format: >ID\\nSEQUENCE or just SEQUENCE):"))
        self._msa_input = QTextEdit()
        self._msa_input.setPlaceholderText(">seq1\nATGGTGCACCTGACT\n>seq2\nATGGTGCATCTGACT\n>seq3\nATGGTGCACCTGGCT")
        self._msa_input.setFont(QFont("JetBrains Mono", 11))
        msa_layout.addWidget(self._msa_input)

        btn_msa = QPushButton("Run MSA (MUSCLE or fallback)")
        btn_msa.setObjectName("primary")
        btn_msa.clicked.connect(self._run_msa)
        msa_layout.addWidget(btn_msa)

        self._msa_progress = QProgressBar()
        self._msa_progress.setRange(0, 0)
        self._msa_progress.hide()
        msa_layout.addWidget(self._msa_progress)

        self._msa_status = QLabel("")
        self._msa_status.setObjectName("subheading")
        msa_layout.addWidget(self._msa_status)

        self._msa_result = QTextEdit()
        self._msa_result.setReadOnly(True)
        self._msa_result.setFont(QFont("JetBrains Mono", 11))
        msa_layout.addWidget(self._msa_result)

        # Identity matrix
        msa_layout.addWidget(QLabel("Pairwise identity matrix (%):"))
        self._identity_table = QTableWidget()
        self._identity_table.setMaximumHeight(150)
        msa_layout.addWidget(self._identity_table)

        tabs.addTab(msa_widget, "Multiple Sequence Alignment")

        layout.addWidget(tabs)

    def _run_pairwise(self) -> None:
        s1 = self._seq1.toPlainText().strip().upper().replace(" ", "").replace("\n", "")
        s2 = self._seq2.toPlainText().strip().upper().replace(" ", "").replace("\n", "")
        if not s1 or not s2:
            QMessageBox.warning(self, "Missing sequences", "Paste both sequences first.")
            return

        mode = self._mode.currentData()
        aligner = Align.PairwiseAligner()
        aligner.mode = mode
        aligner.match_score = 2.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -2.0
        aligner.extend_gap_score = -0.5

        try:
            alignments = list(aligner.align(s1, s2))
            if not alignments:
                self._pair_result.setPlainText("No alignment found.")
                return
            best = alignments[0]
            counts = best.counts()
            aln_len = best.length

            identity = counts.identities / aln_len * 100 if aln_len else 0
            fasta_lines = best.format("fasta").strip().split("\n")
            gapped = [ln for ln in fasta_lines if not ln.startswith(">")]
            a1 = gapped[0] if len(gapped) > 0 else s1
            a2 = gapped[1] if len(gapped) > 1 else s2

            self._stat_score.setText(f"Score: {best.score:.1f}")
            self._stat_identity.setText(f"Identity: {identity:.1f}%")
            self._stat_gaps.setText(f"Gaps: {counts.gaps}")
            self._stat_len.setText(f"Length: {aln_len}")

            # Build visual alignment with match line
            match_line = "".join(
                "|" if a == b else "." if (a != "-" and b != "-") else " "
                for a, b in zip(a1, a2)
            )
            display = f"Seq1  {a1}\n      {match_line}\nSeq2  {a2}\n\nIdentity: {identity:.1f}% | Score: {best.score:.1f} | Gaps: {counts.gaps}"
            self._pair_result.setPlainText(display)

        except Exception as e:
            self._pair_result.setPlainText(f"Alignment error: {e}")

    def _run_msa(self) -> None:
        raw = self._msa_input.toPlainText().strip()
        if not raw:
            return

        # Parse FASTA-like or plain sequences
        seqs = []
        current_id = None
        current_seq = []
        for line in raw.splitlines():
            line = line.strip()
            if line.startswith(">"):
                if current_seq and current_id:
                    seqs.append({"id": current_id, "seq": "".join(current_seq).upper()})
                current_id = line[1:].split()[0]
                current_seq = []
            elif line:
                if current_id is None:
                    current_id = f"seq{len(seqs) + 1}"
                current_seq.append(line)
        if current_seq and current_id:
            seqs.append({"id": current_id, "seq": "".join(current_seq).upper()})

        if len(seqs) < 2:
            QMessageBox.warning(self, "Too few sequences", "Need at least 2 sequences for MSA.")
            return

        self._msa_progress.show()
        self._msa_status.setText(f"Aligning {len(seqs)} sequences…")

        self._msa_worker = Worker(self._do_msa, seqs)
        self._msa_worker.result.connect(self._on_msa_done)
        self._msa_worker.error.connect(lambda e: (
            self._msa_progress.hide(),
            self._msa_status.setText(f"Error: {e}"),
        ))
        self._msa_worker.start()

    def _do_msa(self, seqs: list[dict]) -> dict:
        import subprocess
        import tempfile
        import os
        fasta_in = "".join(f">{s['id']}\n{s['seq']}\n" for s in seqs)
        try:
            with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fin:
                fin.write(fasta_in)
                fin_path = fin.name
            out_path = fin_path + ".aln"
            result = subprocess.run(
                ["muscle", "-align", fin_path, "-output", out_path],
                capture_output=True, timeout=60,
            )
            if result.returncode != 0:
                raise RuntimeError(result.stderr.decode())
            from Bio import SeqIO
            aligned = [(r.id, str(r.seq)) for r in SeqIO.parse(out_path, "fasta")]
            os.unlink(fin_path)
            os.unlink(out_path)
        except (FileNotFoundError, RuntimeError):
            # Fallback: pad to same length
            os.unlink(fin_path) if "fin_path" in dir() and os.path.exists(fin_path) else None
            max_len = max(len(s["seq"]) for s in seqs)
            aligned = [(s["id"], s["seq"] + "-" * (max_len - len(s["seq"]))) for s in seqs]

        # Consensus
        if aligned:
            length = max(len(seq) for _, seq in aligned)
            consensus = ""
            for i in range(length):
                col = [seq[i].upper() for _, seq in aligned if i < len(seq)]
                most = max(set(col), key=col.count) if col else "N"
                consensus += most if col.count(most) > len(col) / 2 else "N"
        else:
            consensus = ""

        # Identity matrix
        n = len(aligned)
        matrix = [[0.0] * n for _ in range(n)]
        for i in range(n):
            matrix[i][i] = 100.0
            for j in range(i + 1, n):
                s1, s2 = aligned[i][1], aligned[j][1]
                same = sum(a == b for a, b in zip(s1, s2) if a != "-" and b != "-")
                total = sum(1 for a, b in zip(s1, s2) if a != "-" or b != "-")
                pct = round(same / total * 100, 1) if total else 0.0
                matrix[i][j] = matrix[j][i] = pct

        return {"aligned": aligned, "consensus": consensus, "matrix": matrix}

    def _on_msa_done(self, result: dict) -> None:
        self._msa_progress.hide()
        aligned = result["aligned"]
        self._msa_status.setText(f"Aligned {len(aligned)} sequences.")

        # Display alignment
        id_width = max(len(i) for i, _ in aligned) + 2
        lines = []
        for seq_id, seq in aligned:
            lines.append(f"{seq_id.ljust(id_width)}{seq}")
        lines.append(f"{'consensus'.ljust(id_width)}{result['consensus']}")
        self._msa_result.setPlainText("\n".join(lines))

        # Identity matrix
        ids = [i for i, _ in aligned]
        n = len(ids)
        mat = result["matrix"]
        self._identity_table.setRowCount(n)
        self._identity_table.setColumnCount(n)
        self._identity_table.setHorizontalHeaderLabels(ids)
        self._identity_table.setVerticalHeaderLabels(ids)
        for i in range(n):
            for j in range(n):
                val = mat[i][j]
                item = QTableWidgetItem(f"{val:.1f}")
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                if val == 100.0:
                    item.setBackground(QColor("#1e3a2e"))
                elif val > 80:
                    item.setBackground(QColor("#1a2d3a"))
                self._identity_table.setItem(i, j, item)
