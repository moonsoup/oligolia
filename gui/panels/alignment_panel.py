"""Sequence alignment panel — pairwise and MSA."""

from __future__ import annotations
import html
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, QPushButton,
    QLabel, QComboBox, QTabWidget, QTableWidgetItem,
    QProgressBar, QMessageBox,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QColor
from Bio import Align
from ..workers import Worker, worker_busy
from ..table_header import fitted_table


#: Columns per block. 60 is the ClustalW/BLAST convention and fits comfortably
#: without horizontal scrolling.
BLOCK_WIDTH = 60


def format_alignment_blocks(a1: str, a2: str, width: int = BLOCK_WIDTH) -> str:
    """Render a pairwise alignment as fixed-width blocks with coordinates.

    The previous rendering put the whole alignment on three very long lines and
    left QTextEdit to word-wrap them. The match line contains spaces (wherever
    either side has a gap) while the sequence lines contain none, so the three
    wrapped at different columns and the bars stopped lining up with the bases
    they describe — correct numbers over a misleading picture (#77).

    Blocks remove the question entirely: every row in a block is exactly `width`
    characters, so the columns cannot drift however the widget is sized.

    Coordinates count bases, not columns, so a gap does not advance them — which
    is what a user needs in order to find a mismatch in the original sequence.
    """
    if len(a1) != len(a2):
        # Shouldn't happen for an aligned pair; don't pretend it lines up.
        raise ValueError(f"aligned sequences differ in length: {len(a1)} vs {len(a2)}")

    label_w = 4
    pos_w = max(len(str(len(a1))), 5)
    out: list[str] = []
    top = bottom = 0  # bases consumed so far, excluding gaps

    for start in range(0, len(a1), width):
        chunk1 = a1[start:start + width]
        chunk2 = a2[start:start + width]
        match = "".join(
            "|" if x == y and x != "-" else "." if (x != "-" and y != "-") else " "
            for x, y in zip(chunk1, chunk2)
        )

        top_start = top + 1
        bottom_start = bottom + 1
        top += len(chunk1) - chunk1.count("-")
        bottom += len(chunk2) - chunk2.count("-")

        # The match row is padded to the same total width as the sequence rows,
        # including the trailing coordinate field. Codex pointed out that leaving
        # it short made "every row in a block is the same width" false, even
        # though the columns themselves lined up.
        out.append(f"{'Seq1':<{label_w}} {top_start:>{pos_w}} {chunk1} {top:>{pos_w}}")
        out.append(f"{'':<{label_w}} {'':>{pos_w}} {match:<{len(chunk1)}} {'':>{pos_w}}")
        out.append(f"{'Seq2':<{label_w}} {bottom_start:>{pos_w}} {chunk2} {bottom:>{pos_w}}")
        out.append("")

    return "\n".join(out).rstrip("\n")


def format_aligner_notice(statuses: list) -> str:
    """One line saying what to install, built from the router's own hint.

    Takes `backend.routers.alignment.AlignerInfo`s. The wording is not restated
    here: the sentence quotes `hint`, the same string `/alignment/multiple`'s 503
    quotes, with the aligner's URL turned into a clickable link (#83).
    """
    from backend.routers.alignment import DEFAULT_ALGORITHM

    primary = next(
        (a for a in statuses if a.name == DEFAULT_ALGORITHM),
        statuses[0] if statuses else None,
    )
    if primary is None:  # no aligners known at all — nothing to advise
        return ""
    url = html.escape(primary.url, quote=True)
    hint = html.escape(primary.hint).replace(
        html.escape(primary.url), f'<a href="{url}">{html.escape(primary.url)}</a>'
    )
    return (
        "No sequence aligner found, so multiple alignment cannot run. "
        f"Install {hint}, then reopen this tab."
    )


class AlignmentPanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self._msa_worker: Worker | None = None
        #: aligner name → its display spelling, as the router reports it. Filled
        #: by `_populate_aligner_choices`; read off the GUI thread is fine, it is
        #: only ever written there.
        self._aligner_labels: dict[str, str] = {}
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
        # Blocks are already fixed-width; wrapping them would undo that (#77).
        self._pair_result.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self._pair_result.setFont(QFont("JetBrains Mono", 11))
        pair_layout.addWidget(self._pair_result)

        tabs.addTab(pair_widget, "Pairwise")

        # ── MSA ───────────────────────────────────────────────────────────────
        msa_widget = QWidget()
        msa_layout = QVBoxLayout(msa_widget)

        msa_layout.addWidget(QLabel("Enter sequences (one per line, format: >ID\\nSEQUENCE or just SEQUENCE):"))

        # Says what to install, before the user pastes anything (#83). Hidden
        # whenever an aligner is present, so the tab looks exactly as it did.
        self._msa_notice = QLabel()
        self._msa_notice.setObjectName("subheading")
        self._msa_notice.setTextFormat(Qt.TextFormat.RichText)
        self._msa_notice.setOpenExternalLinks(True)
        self._msa_notice.setWordWrap(True)
        self._msa_notice.hide()
        msa_layout.addWidget(self._msa_notice)

        # Which aligner actually runs. The router's MSARequest has always taken an
        # `algorithm`; the tab used to pin it to the default, so a machine with
        # only ClustalW installed got enabled controls and then a 503 naming the
        # aligner it does not have. The selection made here is the one passed to
        # `multiple_align` (#83).
        aligner_row = QHBoxLayout()
        self._msa_algorithm = QComboBox()
        self._msa_algorithm.currentIndexChanged.connect(self._sync_msa_button)
        aligner_row.addWidget(QLabel("Aligner:"))
        aligner_row.addWidget(self._msa_algorithm)
        aligner_row.addStretch()
        msa_layout.addLayout(aligner_row)

        self._msa_input = QTextEdit()
        self._msa_input.setPlaceholderText(">seq1\nATGGTGCACCTGACT\n>seq2\nATGGTGCATCTGACT\n>seq3\nATGGTGCACCTGGCT")
        self._msa_input.setFont(QFont("JetBrains Mono", 11))
        msa_layout.addWidget(self._msa_input)

        self._btn_msa = QPushButton("Run MSA (requires MUSCLE)")
        self._btn_msa.setObjectName("primary")
        self._btn_msa.clicked.connect(self._run_msa)
        msa_layout.addWidget(self._btn_msa)

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
        self._identity_table = fitted_table()
        self._identity_table.setMaximumHeight(150)
        msa_layout.addWidget(self._identity_table)

        self._msa_tab_index = tabs.addTab(msa_widget, "Multiple Sequence Alignment")

        layout.addWidget(tabs)
        self._tabs = tabs
        tabs.currentChanged.connect(self._on_tab_changed)
        self._refresh_aligner_availability()

    # ── Aligner availability ─────────────────────────────────────────────────

    def _on_tab_changed(self, index: int) -> None:
        # Re-check on every visit, so someone who installs MUSCLE while the app
        # is open does not have to restart it (#83).
        if index == self._msa_tab_index:
            self._refresh_aligner_availability()

    def showEvent(self, event) -> None:  # noqa: N802 — Qt's spelling
        super().showEvent(event)
        if self._tabs.currentIndex() == self._msa_tab_index:
            self._refresh_aligner_availability()

    def _refresh_aligner_availability(self) -> None:
        """Enable or disable the MSA controls from the router's own detection.

        Imported here rather than at module import for the same reason `_do_msa`
        does it: the panel is built at startup and the router pulls in FastAPI.
        """
        from backend.routers.alignment import (
            DEFAULT_ALGORITHM, aligner_statuses, preferred_algorithm,
        )

        statuses = aligner_statuses()
        ready = any(a.available for a in statuses)
        # Keep the user's own pick across a re-check, but only while it is still
        # runnable — otherwise fall back to the router's preference.
        keep = self._msa_algorithm.currentData()
        if keep is not None and not any(a.name == keep and a.available for a in statuses):
            keep = None
        self._populate_aligner_choices(statuses, keep or preferred_algorithm(statuses))

        self._msa_input.setEnabled(ready)
        self._msa_algorithm.setEnabled(ready)
        self._btn_msa.setEnabled(ready)
        notice = "" if ready else format_aligner_notice(statuses)
        self._msa_notice.setText(notice)
        self._msa_notice.setVisible(not ready)
        # Plain text, because a tooltip does not render the anchor usefully.
        self._btn_msa.setToolTip("" if ready else next(
            (a.hint for a in statuses if a.name == DEFAULT_ALGORITHM), ""
        ))
        self._sync_msa_button()

    def _populate_aligner_choices(self, statuses: list, select: str | None) -> None:
        """Refill the selector: every aligner listed, only installed ones pickable.

        A missing aligner stays visible and greyed rather than disappearing, so
        the list is also the answer to "what could I install?".
        """
        blocked = self._msa_algorithm.blockSignals(True)
        try:
            self._msa_algorithm.clear()
            self._aligner_labels = {a.name: a.label for a in statuses}
            for row, info in enumerate(statuses):
                self._msa_algorithm.addItem(
                    info.label if info.available else f"{info.label} (not installed)",
                    info.name,
                )
                if not info.available:
                    item = self._msa_algorithm.model().item(row)
                    if item is not None:
                        item.setEnabled(False)
                    self._msa_algorithm.setItemData(row, info.hint, Qt.ItemDataRole.ToolTipRole)
            index = self._msa_algorithm.findData(select) if select else -1
            self._msa_algorithm.setCurrentIndex(index)
        finally:
            self._msa_algorithm.blockSignals(blocked)

    def _sync_msa_button(self) -> None:
        """Name the aligner the button will actually run."""
        selected = self._msa_algorithm.currentData()
        self._btn_msa.setText(
            f"Run MSA ({self._aligner_labels.get(selected, selected)})"
            if selected else "Run MSA (requires MUSCLE)"
        )

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
            # Lazily take the first (optimal) alignment. Never list() the result and
            # never len() it: two unrelated 500-nt sequences already have more
            # co-optimal alignments than fit in an int64, which showed up here as a
            # blank "Alignment error: " (#56).
            try:
                best = next(iter(aligner.align(s1, s2)))
            except StopIteration:
                self._pair_result.setPlainText("No alignment found.")
                return
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

            self._pair_result.setPlainText(format_alignment_blocks(a1, a2))

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

        # Read on the GUI thread; the worker must not touch widgets.
        algorithm = self._msa_algorithm.currentData()
        if algorithm is None:
            # Only reachable if detection went stale between the last re-check and
            # the click; re-check rather than shell out to nothing.
            self._refresh_aligner_availability()
            self._msa_status.setText("No sequence aligner is installed.")
            return

        self._msa_progress.show()
        self._msa_status.setText(
            f"Aligning {len(seqs)} sequences with "
            f"{self._aligner_labels.get(algorithm, algorithm)}…"
        )

        # Refuse rather than rebind a live QThread, which Qt aborts on (#54).
        if worker_busy(self, "_msa_worker"):
            self._msa_progress.hide()
            self._msa_status.setText("Already aligning — wait for that run to finish.")
            return
        self._msa_worker = Worker(self._do_msa, seqs, algorithm)
        self._msa_worker.result.connect(self._on_msa_done)
        self._msa_worker.error.connect(lambda e: (
            self._msa_progress.hide(),
            self._msa_status.setText(f"Error: {e}"),
        ))
        self._msa_worker.start()

    def _do_msa(self, seqs: list[dict], algorithm: str | None = None) -> dict:
        """Align via the backend router, in-process, with the selected aligner.

        This used to be a second, independent copy of the MUSCLE call, the
        right-padding fallback, the consensus and the identity matrix — so #58 had
        to be fixed twice and the two could drift apart. It now calls the router,
        which refuses rather than approximating; the Worker's `error` signal puts
        the reason on the status line.

        `algorithm` comes from the selector, which is populated from the same
        detection the router checks before shelling out — so the tab cannot offer
        an aligner this call then 503s on (#83). None means "no preference": ask
        the router which one it would pick, and fall back to its stated default
        so that with nothing installed the refusal still names something.
        """
        from fastapi import HTTPException

        from backend.routers import alignment

        if algorithm is None:
            algorithm = alignment.preferred_algorithm() or alignment.DEFAULT_ALGORITHM

        try:
            res = alignment.multiple_align(
                alignment.MSARequest(sequences=seqs, algorithm=algorithm)
            )
        except HTTPException as e:
            # Worker.error stringifies whatever is raised; HTTPException's own repr
            # buries the message, so surface just the detail.
            raise RuntimeError(e.detail) from None

        return {
            "aligned": [(a["id"], a["aligned_seq"]) for a in res.aligned],
            "consensus": res.consensus,
            "matrix": res.identity_matrix,
        }

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
