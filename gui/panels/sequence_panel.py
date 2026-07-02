"""Sequence viewer and editor panel."""

from __future__ import annotations
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QListWidget, QListWidgetItem,
    QTextEdit, QLabel, QPushButton, QComboBox, QLineEdit, QGroupBox,
    QFormLayout, QSpinBox, QMessageBox, QFileDialog, QApplication,
    QToolBar, QCheckBox, QTableWidget, QTableWidgetItem, QHeaderView,
    QDialog, QDialogButtonBox,
)
from PyQt6.QtCore import Qt, QSettings, pyqtSignal
from datetime import datetime, timezone
from PyQt6.QtGui import (
    QFont, QColor, QTextCharFormat, QSyntaxHighlighter, QTextCursor,
    QKeySequence, QShortcut, QPainter, QPen,
)

from Bio.Seq import Seq
from backend.models.sequence import Sequence, MoleculeType, Annotation
from backend.formats import read_fasta, read_fastq, read_genbank, read_snapgene, VENDORS
from gui.history import UndoStack
from gui.panels.feature_colors import feature_color_map
from gui.panels.plasmid_map import PlasmidMapWidget
import re

# Edit operations that mutate the active sequence in place (and are undoable).
# Molecule-type-changing ops (translate/transcribe/back_transcribe) stay
# non-destructive and render into the result pane instead.
IN_PLACE_OPS = {"insert", "delete", "replace", "reverse_complement", "complement"}


class DNAHighlighter(QSyntaxHighlighter):
    """
    Color-codes ALL IUPAC nucleotide symbols (NCBI/IUPAC 1985 standard).

    Definite bases:
      A=green  T/U=red  G=gold  C=blue

    Ambiguity codes (two-base):
      R(A|G)=lime  Y(C|T)=salmon  M(A|C)=teal  K(G|T)=orange
      W(A|T)=khaki  S(G|C)=purple

    Ambiguity codes (three-base):
      B(not A)=pink  D(not C)=sienna  H(not G)=peru  V(not T)=orchid

    Universal:
      N=gray  -(gap)=dark-gray  .(gap)=dark-gray
    """
    COLORS: dict[str, str] = {
        # Definite bases
        "A": "#4ade80",   # green
        "T": "#f87171",   # red
        "U": "#fb7185",   # rose (RNA uracil)
        "G": "#fbbf24",   # amber/gold
        "C": "#60a5fa",   # blue
        # Two-base ambiguities
        "R": "#a3e635",   # lime       A or G (puRine)
        "Y": "#fca5a5",   # light-red  C or T (pYrimidine)
        "M": "#2dd4bf",   # teal       A or C (aMino)
        "K": "#fb923c",   # orange     G or T (Keto)
        "W": "#d9f99d",   # yellow-grn A or T (Weak)
        "S": "#a78bfa",   # violet     G or C (Strong)
        # Three-base ambiguities
        "B": "#f9a8d4",   # pink       C, G or T (not A)
        "D": "#fcd34d",   # yellow     A, G or T (not C)
        "H": "#fdba74",   # peach      A, C or T (not G)
        "V": "#c4b5fd",   # lavender   A, C or G (not T/U)
        # Universal
        "N": "#94a3b8",   # slate-gray (aNy)
        "-": "#475569",   # dark-gray  (gap)
        ".": "#334155",   # darker-gray (alignment gap)
    }

    def highlightBlock(self, text: str) -> None:
        for i, ch in enumerate(text.upper()):
            color = self.COLORS.get(ch)
            if color:
                fmt = QTextCharFormat()
                fmt.setForeground(QColor(color))
                self.setFormat(i, 1, fmt)


class ProteinHighlighter(QSyntaxHighlighter):
    """
    Color-codes ALL IUPAC amino acid symbols by biochemical property.

    Standard 20 + selenocysteine (U) + pyrrolysine (O) + ambiguity (B/Z/J/X) + stop (*).

    Color scheme (ClustalX-inspired with biochemical grouping):
      Hydrophobic nonpolar  (A V I L M F W P): amber/yellow
      Aromatic              (F W Y H):          orange
      Positively charged    (K R H):            blue
      Negatively charged    (D E):              red
      Polar uncharged       (S T N Q):          green
      Special/sulfur        (C U):              yellow-green (Cys/Sec have SH)
      Glycine               (G):                white (flexible)
      Proline               (P):                purple (rigid)
      Pyrrolysine           (O):                teal
      Ambiguous             (B Z J X):          slate-gray
      Stop codon            (*):                bright red background
      Gap                   (- .):              dark-gray
    """
    # (foreground_color, bold)
    COLORS: dict[str, tuple[str, bool]] = {
        # Hydrophobic nonpolar
        "A": ("#fbbf24", False),   # amber
        "V": ("#f59e0b", False),
        "I": ("#d97706", False),
        "L": ("#b45309", False),
        "M": ("#92400e", False),
        # Aromatic
        "F": ("#fb923c", True),    # orange
        "W": ("#ea580c", True),
        "Y": ("#c2410c", True),
        # Positively charged
        "K": ("#60a5fa", True),    # blue
        "R": ("#3b82f6", True),
        "H": ("#93c5fd", False),   # light-blue (partial positive at pH 7)
        # Negatively charged
        "D": ("#f87171", True),    # red
        "E": ("#dc2626", True),
        # Polar uncharged
        "S": ("#4ade80", False),   # green
        "T": ("#22c55e", False),
        "N": ("#16a34a", False),
        "Q": ("#15803d", False),
        # Special/sulfur-containing
        "C": ("#a3e635", True),    # yellow-green (Cys – forms disulfide)
        "U": ("#84cc16", True),    # Selenocysteine
        # Glycine (flexible – no side chain)
        "G": ("#e2e8f0", False),   # near-white
        # Proline (rigid – cyclized backbone)
        "P": ("#a78bfa", True),    # violet
        # Pyrrolysine (22nd amino acid)
        "O": ("#2dd4bf", False),   # teal
        # Ambiguous
        "B": ("#94a3b8", False),   # Asp or Asn
        "Z": ("#94a3b8", False),   # Glu or Gln
        "J": ("#94a3b8", False),   # Leu or Ile
        "X": ("#64748b", False),   # any
        # Gap
        "-": ("#475569", False),
        ".": ("#334155", False),
    }
    STOP_COLOR = QColor("#7f1d1d")  # dark-red background for stop codon *

    def highlightBlock(self, text: str) -> None:
        for i, ch in enumerate(text.upper()):
            if ch == "*":
                fmt = QTextCharFormat()
                fmt.setBackground(self.STOP_COLOR)
                fmt.setForeground(QColor("#fca5a5"))
                fmt.setFontWeight(700)
                self.setFormat(i, 1, fmt)
            elif ch in self.COLORS:
                color, bold = self.COLORS[ch]
                fmt = QTextCharFormat()
                fmt.setForeground(QColor(color))
                if bold:
                    fmt.setFontWeight(700)
                self.setFormat(i, 1, fmt)


class ReadingFrameHighlighter(QSyntaxHighlighter):
    """
    Background-shades codons (3-base groups) of one reading frame, layered on
    top of the per-base foreground highlighter on the same document.

    Forward frames (+1/+2/+3) shade directly from the displayed sequence.
    Reverse frames (-1/-2/-3) are computed against the reverse complement and
    mapped back onto the forward-oriented display positions.
    """
    BAND_COLORS = ["#1e3a5f", "#3f2d1e", "#1e3f2d"]  # cycles every codon

    def __init__(self, document) -> None:
        super().__init__(document)
        self._frame: int | None = None
        self._seq: str = ""

    def setFrame(self, frame: int | None, seq: str) -> None:
        self._frame = frame
        self._seq = seq
        self.rehighlight()

    def highlightBlock(self, text: str) -> None:
        if not self._frame or not self._seq or len(text) != len(self._seq):
            return
        n = len(self._seq)
        offset = abs(self._frame) - 1
        if self._frame > 0:
            for i in range(offset, n):
                fmt = QTextCharFormat()
                fmt.setBackground(QColor(self.BAND_COLORS[(i - offset) // 3 % 3]))
                self.setFormat(i, 1, fmt)
        else:
            for rc_i in range(offset, n):
                orig_i = n - 1 - rc_i
                fmt = QTextCharFormat()
                fmt.setBackground(QColor(self.BAND_COLORS[(rc_i - offset) // 3 % 3]))
                self.setFormat(orig_i, 1, fmt)


class FeatureHighlighter(QSyntaxHighlighter):
    """Background-shades annotated gene features (CDS, exon, promoter, …) by type."""

    def __init__(self, document) -> None:
        super().__init__(document)
        self._annotations: list[Annotation] = []
        self._color_map: dict[str, QColor] = {}

    def setAnnotations(self, annotations: list[Annotation]) -> None:
        self._annotations = annotations or []
        self._color_map = feature_color_map(self._annotations)
        self.rehighlight()

    @property
    def color_map(self) -> dict[str, QColor]:
        return self._color_map

    def highlightBlock(self, text: str) -> None:
        if not self._annotations:
            return
        for ann in self._annotations:
            start, end = max(0, ann.start), min(len(text), ann.end)
            if start >= end:
                continue
            fmt = QTextCharFormat()
            fmt.setBackground(self._color_map[ann.feature_type])
            # Auto-detected features (issue #44) get a dashed underline so they
            # read as distinct from user/GenBank-provided annotations.
            if ann.qualifiers.get("auto_detected") == "true":
                fmt.setUnderlineStyle(QTextCharFormat.UnderlineStyle.DashUnderline)
                fmt.setUnderlineColor(QColor("#e2e8f0"))
            self.setFormat(start, end - start, fmt)


class FeatureMinimap(QWidget):
    """Slim vertical map of feature ticks along the full sequence length (#33).

    One colored tick per annotation at ``start / length`` down the widget,
    independent of scroll position, so off-screen features stay visible.
    Clicking a tick emits the feature's span for the panel to jump to.
    """

    feature_clicked = pyqtSignal(int, int)  # (start, end)
    _MARGIN = 6

    def __init__(self) -> None:
        super().__init__()
        self._annotations: list[Annotation] = []
        self._length = 1
        self._color_map: dict[str, QColor] = {}
        self.setFixedWidth(22)
        self.setToolTip("Feature map — click a tick to jump to that feature")

    def set_features(self, annotations: list[Annotation], length: int) -> None:
        self._annotations = annotations or []
        self._length = max(length, 1)
        self._color_map = feature_color_map(self._annotations)
        self.update()

    def _y_for(self, position: int) -> int:
        usable = max(self.height() - 2 * self._MARGIN, 1)
        frac = min(max(position / self._length, 0.0), 1.0)
        return self._MARGIN + int(frac * usable)

    def paintEvent(self, event) -> None:  # noqa: N802 (Qt override)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        cx = self.width() // 2
        painter.setPen(QPen(QColor("#334155"), 2))
        painter.drawLine(cx, self._MARGIN, cx, self.height() - self._MARGIN)
        for ann in self._annotations:
            painter.setPen(QPen(self._color_map.get(ann.feature_type, QColor("#94a3b8")), 3))
            y = self._y_for(ann.start)
            painter.drawLine(cx - 7, y, cx + 7, y)
        painter.end()

    def mousePressEvent(self, event) -> None:  # noqa: N802 (Qt override)
        if not self._annotations:
            return
        y = event.position().y()
        nearest = min(self._annotations, key=lambda a: abs(self._y_for(a.start) - y))
        if abs(self._y_for(nearest.start) - y) <= 10:  # click near a tick
            self.feature_clicked.emit(nearest.start, nearest.end)


# Bump this when TERMS.md changes materially so users re-accept (issue #28).
# Keep in sync with the "terms-version" marker in TERMS.md.
_TERMS_VERSION = "1"
_TERMS_URL = "https://github.com/moonsoup/oligolia/blob/main/TERMS.md"


class TermsDialog(QDialog):
    """One-time Terms-of-Use acceptance for synthesis export (issue #28)."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Export for Synthesis — Terms of Use")
        self.setMinimumWidth(460)
        layout = QVBoxLayout(self)

        summary = QLabel(
            "<b>Before exporting a synthesis order file, please confirm:</b>"
            "<ul>"
            "<li>This <b>generates a file only</b> — Oligolia does <b>not</b> submit, "
            "transmit, or place any order with any vendor.</li>"
            "<li>You submit the order yourself by uploading the file on the "
            "vendor's own portal, under their terms and pricing.</li>"
            "<li>Oligolia does <b>not</b> verify sequence correctness, "
            "manufacturability, legality, or safety — <b>review every sequence "
            "before ordering</b> and comply with biosecurity screening and law.</li>"
            "</ul>"
            f'Full terms: <a href="{_TERMS_URL}">TERMS.md</a>'
        )
        summary.setWordWrap(True)
        summary.setTextFormat(Qt.TextFormat.RichText)
        summary.setOpenExternalLinks(True)
        layout.addWidget(summary)

        self._accept_check = QCheckBox("I have read and accept the Terms of Use")
        self._accept_check.toggled.connect(self._on_toggle)
        layout.addWidget(self._accept_check)

        self._buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel | QDialogButtonBox.StandardButton.Ok)
        self._buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Accept and Export")
        self._buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(False)
        self._buttons.accepted.connect(self.accept)
        self._buttons.rejected.connect(self.reject)
        layout.addWidget(self._buttons)

    def _on_toggle(self, checked: bool) -> None:
        # Explicit affirmative action required — OK stays disabled until checked.
        self._buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(checked)


class SequencePanel(QWidget):
    sequence_selected = pyqtSignal(object)  # emits Sequence

    def __init__(self) -> None:
        super().__init__()
        self._sequences: dict[str, Sequence] = {}
        self._active: Sequence | None = None
        self._history: dict[str, UndoStack] = {}  # per-sequence undo/redo
        self._build_ui()
        self._install_shortcuts()

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

        # ── Viewport options bar: selection + view/format toggles ─────────────
        options_bar = QToolBar()
        options_bar.setMovable(False)

        # Topology toggle — always visible, editable per selected sequence.
        self._circular_check = QCheckBox("Circular")
        self._circular_check.setToolTip(
            "Circular (plasmid) vs linear topology. Affects origin-spanning "
            "restriction/digest analysis; preserved across edits and GenBank export.")
        self._circular_check.toggled.connect(self._on_circular_toggled)
        options_bar.addWidget(self._circular_check)
        options_bar.addSeparator()

        options_bar.addWidget(QLabel(" Select: "))
        self._sel_start = QSpinBox(); self._sel_start.setMaximum(999_999_999)
        self._sel_end = QSpinBox(); self._sel_end.setMaximum(999_999_999)
        options_bar.addWidget(self._sel_start)
        options_bar.addWidget(QLabel("–"))
        options_bar.addWidget(self._sel_end)
        btn_select = QPushButton("Select")
        btn_select.clicked.connect(self._select_range)
        options_bar.addWidget(btn_select)

        options_bar.addSeparator()

        options_bar.addWidget(QLabel(" Reading frame: "))
        self._frame_combo = QComboBox()
        self._frame_combo.addItem("Off", None)
        for f in (1, 2, 3, -1, -2, -3):
            self._frame_combo.addItem(f"{f:+d}", f)
        self._frame_combo.currentIndexChanged.connect(self._on_frame_changed)
        options_bar.addWidget(self._frame_combo)

        options_bar.addSeparator()

        self._feature_check = QCheckBox("Highlight features")
        self._feature_check.toggled.connect(self._on_features_toggled)
        options_bar.addWidget(self._feature_check)

        # Plasmid-map option: unique cutters only (shown only for circular seqs).
        self._unique_cutters_check = QCheckBox("Unique cutters")
        self._unique_cutters_check.setChecked(True)
        self._unique_cutters_check.setToolTip("Plasmid map: show only single-cutter enzymes")
        self._unique_cutters_check.toggled.connect(
            lambda c: self._plasmid_map.set_unique_cutters_only(c))
        self._unique_cutters_check.hide()
        options_bar.addWidget(self._unique_cutters_check)

        self._feature_legend = QLabel("")
        self._feature_legend.setTextFormat(Qt.TextFormat.RichText)
        options_bar.addWidget(self._feature_legend)

        right_layout.addWidget(options_bar)

        # Sequence display
        self._seq_display = QTextEdit()
        self._seq_display.setReadOnly(True)
        self._seq_display.setFont(QFont("JetBrains Mono,Fira Code,Courier New", 11))
        # Per-base highlighter is set dynamically per molecule type (see _apply_highlighter).
        # Frame/feature highlighters are layered on top and stay attached; they no-op when off.
        self._dna_highlighter = DNAHighlighter(self._seq_display.document())
        self._prot_highlighter: ProteinHighlighter | None = None
        self._frame_highlighter = ReadingFrameHighlighter(self._seq_display.document())
        self._feature_highlighter = FeatureHighlighter(self._seq_display.document())
        # Vertical feature minimap beside the sequence view (issue #33).
        self._minimap = FeatureMinimap()
        self._minimap.feature_clicked.connect(self._jump_to_span)
        linear_view = QWidget()
        lv_layout = QHBoxLayout(linear_view)
        lv_layout.setContentsMargins(0, 0, 0, 0)
        lv_layout.addWidget(self._minimap)
        lv_layout.addWidget(self._seq_display)

        # Circular sequences get a plasmid map alongside the linear view
        # (dual-viewer pattern); linear sequences hide it entirely (#27).
        self._plasmid_map = PlasmidMapWidget()
        self._plasmid_map.hide()
        # Bidirectional sync (#26): map arc click -> select in sequence view;
        # sequence-view selection -> glow the overlapping map arc(s).
        self._plasmid_map.feature_clicked.connect(self._jump_to_span)
        self._seq_display.selectionChanged.connect(self._on_seq_selection_changed)

        viewer_splitter = QSplitter(Qt.Orientation.Horizontal)
        viewer_splitter.addWidget(linear_view)
        viewer_splitter.addWidget(self._plasmid_map)
        viewer_splitter.setSizes([520, 360])
        right_layout.addWidget(viewer_splitter)

        # Feature table (issue #31) — browse the selected sequence's annotations.
        feat_grp = QGroupBox("Features")
        feat_layout = QVBoxLayout(feat_grp)
        self._feature_empty = QLabel("No features annotated")
        self._feature_empty.setObjectName("subheading")
        feat_layout.addWidget(self._feature_empty)
        self._feature_table = QTableWidget()
        self._feature_table.setColumnCount(5)
        self._feature_table.setHorizontalHeaderLabels(
            ["Feature Type", "Strand", "Start", "End", "Qualifiers"])
        self._feature_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)
        self._feature_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._feature_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self._feature_table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self._feature_table.itemSelectionChanged.connect(self._on_feature_row_selected)
        self._feature_table.setSortingEnabled(True)
        self._feature_table.setMaximumHeight(160)
        self._feature_table.hide()
        feat_layout.addWidget(self._feature_table)
        right_layout.addWidget(feat_grp)

        # Result display — separate document with its own highlighter
        self._result_display = QTextEdit()
        self._result_display.setReadOnly(True)
        self._result_display.setPlaceholderText("Operation results appear here…")
        self._result_display.setMaximumHeight(130)
        self._result_display.setFont(QFont("JetBrains Mono,Fira Code,Courier New", 11))
        self._result_dna_hl = DNAHighlighter(self._result_display.document())
        self._result_prot_hl: ProteinHighlighter | None = None

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

        self._btn_undo = QPushButton("↶ Undo")
        self._btn_undo.setToolTip("Undo last edit (Ctrl+Z)")
        self._btn_undo.clicked.connect(self._undo)
        self._btn_undo.setEnabled(False)
        op_row.addWidget(self._btn_undo)

        self._btn_redo = QPushButton("↷ Redo")
        self._btn_redo.setToolTip("Redo (Ctrl+Y / Ctrl+Shift+Z)")
        self._btn_redo.clicked.connect(self._redo)
        self._btn_redo.setEnabled(False)
        op_row.addWidget(self._btn_redo)

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

        # Synthesis vendor export — generates a vendor-formatted order file
        # (uses the Select range above, full sequence if none chosen). No
        # order is submitted; the user uploads the file on the vendor's own portal.
        synth_grp = QGroupBox("Export for Synthesis")
        synth_layout = QHBoxLayout(synth_grp)
        self._vendor_combo = QComboBox()
        for key, profile in VENDORS.items():
            label = profile.display_name + ("" if profile.verified else " — unverified format")
            self._vendor_combo.addItem(label, key)
        synth_layout.addWidget(self._vendor_combo, 1)
        btn_export_synth = QPushButton("Export…")
        btn_export_synth.clicked.connect(self._export_synthesis)
        synth_layout.addWidget(btn_export_synth)
        right_layout.addWidget(synth_grp)

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

    def _apply_highlighter(self, display: QTextEdit, is_protein: bool) -> None:
        """Swap the syntax highlighter on a display widget based on molecule type."""
        if is_protein:
            # Replace DNA highlighter with protein highlighter
            self._dna_highlighter.setDocument(None)
            if self._prot_highlighter is None:
                self._prot_highlighter = ProteinHighlighter(display.document())
            else:
                self._prot_highlighter.setDocument(display.document())
        else:
            if self._prot_highlighter:
                self._prot_highlighter.setDocument(None)
            self._dna_highlighter.setDocument(display.document())

    def _select_span(self, start: int, end: int) -> None:
        """Select [start, end) in the sequence view and scroll it into view."""
        if start > end:
            start, end = end, start
        cursor = self._seq_display.textCursor()
        cursor.setPosition(start)
        cursor.setPosition(end, QTextCursor.MoveMode.KeepAnchor)
        self._seq_display.setTextCursor(cursor)
        self._seq_display.ensureCursorVisible()

    def _select_range(self) -> None:
        if not self._active:
            return
        self._select_span(self._sel_start.value(), self._sel_end.value())
        self._seq_display.setFocus()

    def _on_seq_selection_changed(self) -> None:
        """Glow the map arc(s) overlapping the current sequence-view selection."""
        cursor = self._seq_display.textCursor()
        if cursor.hasSelection():
            self._plasmid_map.set_highlight((cursor.selectionStart(), cursor.selectionEnd()))
        else:
            self._plasmid_map.set_highlight(None)

    def _jump_to_span(self, start: int, end: int) -> None:
        """Select+scroll to [start, end) and mirror it in the Select spinboxes."""
        if not self._active:
            return
        self._select_span(int(start), int(end))
        self._sel_start.setValue(int(start))
        self._sel_end.setValue(int(end))

    def _on_feature_row_selected(self) -> None:
        """Clicking a feature row selects/scrolls to its span (issue #32)."""
        rows = self._feature_table.selectionModel().selectedRows()
        if not rows or not self._active:
            return
        r = rows[0].row()
        start = self._feature_table.item(r, 2).data(Qt.ItemDataRole.DisplayRole)
        end = self._feature_table.item(r, 3).data(Qt.ItemDataRole.DisplayRole)
        self._jump_to_span(int(start), int(end))

    def _on_frame_changed(self) -> None:
        frame = self._frame_combo.currentData()
        seq = self._active.seq if (self._active and frame) else ""
        self._frame_highlighter.setFrame(frame, seq)

    def _on_features_toggled(self, checked: bool) -> None:
        annotations = self._active.annotations if (self._active and checked) else []
        self._feature_highlighter.setAnnotations(annotations)
        self._update_feature_legend()

    def _update_feature_legend(self) -> None:
        if not self._feature_check.isChecked() or not self._feature_highlighter.color_map:
            self._feature_legend.setText("")
            return
        swatches = " ".join(
            f"<span style='background-color:{color.name()};'>&nbsp;&nbsp;</span> {ftype}"
            for ftype, color in self._feature_highlighter.color_map.items()
        )
        self._feature_legend.setText(swatches)

    def _update_op_fields(self) -> None:
        op = self._op_combo.currentData()
        self._pos_widget.setVisible(op in ("delete", "replace"))
        self._insert_widget.setVisible(op in ("insert", "replace"))

    def add_sequence(self, seq: Sequence) -> None:
        self._auto_annotate(seq)
        self._sequences[seq.id] = seq
        item = QListWidgetItem(f"{'🔵' if seq.molecule_type == MoleculeType.DNA else '🟡' if seq.molecule_type == MoleculeType.RNA else '🟣'} {seq.id}")
        item.setData(Qt.ItemDataRole.UserRole, seq.id)
        item.setToolTip(f"{seq.length:,} bp · {seq.molecule_type.value}")
        self._list.addItem(item)
        self._list.setCurrentItem(item)

    # Skip the scan on very large sequences (genomic, not plasmid-scale) to
    # keep loading responsive — the reference parts target plasmids.
    _AUTO_ANNOTATE_MAX_BP = 50_000

    def _auto_annotate(self, seq: Sequence) -> None:
        """Fill in common-part annotations on load when a DNA sequence has none.

        Only runs when there are no existing annotations (never overwrites
        real ones, e.g. from a GenBank file) and never blocks loading on error.
        """
        if seq.annotations or seq.molecule_type != MoleculeType.DNA:
            return
        if len(seq.seq) > self._AUTO_ANNOTATE_MAX_BP:
            return
        try:
            from backend.formats import auto_annotate
            hits = auto_annotate(seq.seq)
        except Exception:
            return  # annotation is best-effort; loading must not fail
        if hits:
            seq.annotations = hits

    def _on_seq_selected(self, item: QListWidgetItem | None, _: object = None) -> None:
        if not item:
            return
        seq_id = item.data(Qt.ItemDataRole.UserRole)
        seq = self._sequences.get(seq_id)
        if not seq:
            return
        self._active = seq
        self._render_active()
        self._update_undo_buttons()
        self.sequence_selected.emit(seq)

    def _render_active(self) -> None:
        """Refresh the viewer widgets from ``self._active`` (no signal emit)."""
        seq = self._active
        if not seq:
            return
        is_protein = seq.molecule_type == MoleculeType.PROTEIN
        gc_str = f"GC: {self._gc(seq.seq):.1f}%" if not is_protein else f"AA: {seq.length}"
        topo = "" if is_protein else f"  ·  {'circular' if seq.is_circular else 'linear'}"
        self._info_label.setText(
            f"{seq.name or seq.id}  ·  {seq.molecule_type.value}  ·  "
            f"{seq.length:,} {'aa' if is_protein else 'bp'}  ·  {gc_str}{topo}"
        )
        # Sync the topology toggle without triggering write-back; disabled for protein.
        self._circular_check.blockSignals(True)
        self._circular_check.setChecked(seq.is_circular)
        self._circular_check.setEnabled(not is_protein)
        self._circular_check.blockSignals(False)
        self._apply_highlighter(self._seq_display, is_protein)
        self._seq_display.setPlainText(seq.seq)
        self._sel_start.setMaximum(seq.length)
        self._sel_end.setMaximum(seq.length)
        self._frame_combo.setEnabled(not is_protein)
        if is_protein:
            self._frame_combo.setCurrentIndex(0)  # "Off" — frames are meaningless for protein
        self._on_frame_changed()
        self._feature_highlighter.setAnnotations(seq.annotations if self._feature_check.isChecked() else [])
        self._update_feature_legend()
        self._populate_feature_table()
        self._minimap.set_features(seq.annotations, seq.length)

        # Plasmid map: shown alongside the linear view for circular DNA only.
        show_map = seq.is_circular and not is_protein
        if show_map:
            self._plasmid_map.set_unique_cutters_only(self._unique_cutters_check.isChecked())
            self._plasmid_map.set_sequence(seq)
        self._plasmid_map.setVisible(show_map)
        self._unique_cutters_check.setVisible(show_map)

    def _populate_feature_table(self) -> None:
        """Fill the feature table from the active sequence's annotations."""
        anns = self._active.annotations if self._active else []
        tbl = self._feature_table
        tbl.blockSignals(True)  # avoid spurious selection signals while rebuilding
        tbl.setSortingEnabled(False)
        tbl.setRowCount(0)
        if not anns:
            tbl.hide()
            self._feature_empty.show()
            tbl.setSortingEnabled(True)
            tbl.blockSignals(False)
            return
        self._feature_empty.hide()
        tbl.show()
        cmap = feature_color_map(anns)  # same colors as FeatureHighlighter
        tbl.setRowCount(len(anns))
        for i, a in enumerate(anns):
            label = (a.qualifiers.get("label") or a.qualifiers.get("gene")
                     or a.qualifiers.get("product") or a.qualifiers.get("note") or "")
            summary = str(label)
            if a.qualifiers.get("auto_detected") == "true":
                summary = (summary + "  ·  auto-detected").strip()

            type_item = QTableWidgetItem(a.feature_type)
            type_item.setBackground(cmap[a.feature_type])
            strand_item = QTableWidgetItem(a.strand.value)
            # Numeric sort for coordinates: store ints in the display role.
            start_item = QTableWidgetItem()
            start_item.setData(Qt.ItemDataRole.DisplayRole, a.start)
            end_item = QTableWidgetItem()
            end_item.setData(Qt.ItemDataRole.DisplayRole, a.end)
            for col, item in enumerate(
                    [type_item, strand_item, start_item, end_item, QTableWidgetItem(summary)]):
                tbl.setItem(i, col, item)
        tbl.setSortingEnabled(True)
        tbl.blockSignals(False)

    def _on_circular_toggled(self, checked: bool) -> None:
        """Persist the topology change onto the active sequence + list tooltip."""
        if not self._active:
            return
        self._active.is_circular = checked
        item = self._list.currentItem()
        if item:
            topo = "circular" if checked else "linear"
            item.setToolTip(f"{self._active.length:,} bp · {self._active.molecule_type.value} · {topo}")
        self._render_active()

    # ── Undo/redo ────────────────────────────────────────────────────────────
    def _install_shortcuts(self) -> None:
        undo_sc = QShortcut(QKeySequence(QKeySequence.StandardKey.Undo), self)
        undo_sc.activated.connect(self._undo)
        redo_sc = QShortcut(QKeySequence(QKeySequence.StandardKey.Redo), self)
        redo_sc.activated.connect(self._redo)
        # Ctrl+Shift+Z as an explicit redo alias (common on macOS/Linux).
        redo_alt = QShortcut(QKeySequence("Ctrl+Shift+Z"), self)
        redo_alt.activated.connect(self._redo)

    def _history_for(self, seq_id: str) -> UndoStack:
        stack = self._history.get(seq_id)
        if stack is None:
            stack = UndoStack()
            self._history[seq_id] = stack
        return stack

    def _set_active_seq(self, new_seq: str) -> None:
        """Replace the active sequence's bases (no history push) and re-render."""
        seq = self._active
        seq.seq = new_seq
        seq.length = len(new_seq)
        self._sequences[seq.id] = seq
        item = self._list.currentItem()
        if item:
            item.setToolTip(f"{seq.length:,} bp · {seq.molecule_type.value}")
        self._render_active()

    def _commit_edit(self, new_seq: str, msg: str) -> None:
        """Apply an in-place edit with an undo checkpoint of the prior state."""
        self._history_for(self._active.id).push(self._active.seq)
        self._set_active_seq(new_seq)
        self._set_result_highlighter(self._active.molecule_type == MoleculeType.PROTEIN)
        self._result_display.setPlainText(f"// {msg}\n{new_seq}")
        self._last_result = new_seq
        self._last_result_id = self._active.id
        self._update_undo_buttons()

    def _undo(self) -> None:
        if not self._active:
            return
        prev = self._history_for(self._active.id).undo(self._active.seq)
        if prev is None:
            return
        self._set_active_seq(prev)
        self._update_undo_buttons()

    def _redo(self) -> None:
        if not self._active:
            return
        nxt = self._history_for(self._active.id).redo(self._active.seq)
        if nxt is None:
            return
        self._set_active_seq(nxt)
        self._update_undo_buttons()

    def _update_undo_buttons(self) -> None:
        stack = self._history.get(self._active.id) if self._active else None
        self._btn_undo.setEnabled(bool(stack and stack.can_undo()))
        self._btn_redo.setEnabled(bool(stack and stack.can_redo()))

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

            # In-place edits mutate the active sequence and are undoable;
            # type-changing ops stay non-destructive in the result pane.
            if op in IN_PLACE_OPS:
                self._commit_edit(result, msg)
                return

            # Switch result-display highlighter based on what the operation produces
            self._set_result_highlighter(op == "translate")
            self._result_display.setPlainText(f"// {msg}\n{result}")
            self._last_result = result
            self._last_result_id = f"{self._active.id}_{op}"

        except Exception as e:
            QMessageBox.critical(self, "Operation failed", str(e))

    def _set_result_highlighter(self, is_protein: bool) -> None:
        """Attach the DNA or protein highlighter to the result-pane document."""
        if is_protein:
            self._result_dna_hl.setDocument(None)
            if self._result_prot_hl is None:
                self._result_prot_hl = ProteinHighlighter(self._result_display.document())
            else:
                self._result_prot_hl.setDocument(self._result_display.document())
        else:
            if self._result_prot_hl:
                self._result_prot_hl.setDocument(None)
            self._result_dna_hl.setDocument(self._result_display.document())

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

    def _ensure_terms_accepted(self) -> bool:
        """Gate synthesis export behind one-time Terms acceptance (issue #28).

        Returns True if the current terms version has already been accepted or
        the user accepts now; False if they decline. Re-prompts only when
        _TERMS_VERSION changes.
        """
        settings = QSettings("Oligolia", "Oligolia")
        if settings.value("synthesis_terms_version", "") == _TERMS_VERSION:
            return True
        if not self._show_terms_dialog():
            return False
        settings.setValue("synthesis_terms_version", _TERMS_VERSION)
        settings.setValue("synthesis_terms_accepted_at", datetime.now(timezone.utc).isoformat())
        return True

    def _show_terms_dialog(self) -> bool:
        return TermsDialog(self).exec() == QDialog.DialogCode.Accepted

    def _export_synthesis(self) -> None:
        # Blocked pending human + attorney review of TERMS.md — #28 shipped a
        # draft terms document that was never reviewed; nothing ships live
        # here until that review completes (see TERMS.md and #29). The real
        # export logic this replaced is in git history on this commit's
        # parent; restore it once terms are actually approved rather than
        # re-deriving it from scratch.
        QMessageBox.information(
            self, "Export for Synthesis — temporarily unavailable",
            "This feature is temporarily disabled while its Terms of Use "
            "undergo human and legal review (see TERMS.md and issue #29). "
            "It will be re-enabled once that review is complete."
        )

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
            "All bioinformatics (*.fasta *.fa *.fna *.faa *.fastq *.fq *.gb *.gbk *.embl *.dna);;"
            "FASTA (*.fasta *.fa *.fna *.faa);;"
            "FASTQ (*.fastq *.fq);;"
            "GenBank (*.gb *.gbk);;"
            "SnapGene (*.dna);;"
            "EMBL (*.embl);;All files (*)"
        )
        if not path:
            return
        ext = path.rsplit(".", 1)[-1].lower()
        try:
            readers = {"fasta": read_fasta, "fa": read_fasta, "fna": read_fasta,
                       "faa": read_fasta, "fastq": read_fastq, "fq": read_fastq,
                       "gb": read_genbank, "gbk": read_genbank, "embl": read_genbank,
                       "dna": read_snapgene}
            reader = readers.get(ext, read_fasta)
            # SnapGene .dna is binary; the rest are text.
            mode = "rb" if ext == "dna" else "r"
            with open(path, mode) as f:
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
            self._history.pop(seq_id, None)
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
