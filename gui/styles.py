"""Dark theme stylesheet for Oligolia PyQt6 GUI."""

DARK_STYLESHEET = """
QMainWindow, QDialog {
    background-color: #0a0f1a;
    color: #e2e8f0;
}

QWidget {
    background-color: #0f172a;
    color: #e2e8f0;
    font-family: "Inter", sans-serif;
    font-size: 13px;
}

QTabWidget::pane {
    border: 1px solid #1e293b;
    background-color: #0f172a;
}

QTabBar::tab {
    background-color: #1e293b;
    color: #94a3b8;
    padding: 8px 16px;
    border: none;
    min-width: 90px;
}

QTabBar::tab:selected {
    background-color: #0f172a;
    color: #4ade80;
    border-bottom: 2px solid #4ade80;
}

QTabBar::tab:hover:!selected {
    background-color: #1e293b;
    color: #cbd5e1;
}

QTextEdit, QPlainTextEdit {
    background-color: #020617;
    color: #4ade80;
    border: 1px solid #1e293b;
    border-radius: 4px;
    font-family: "JetBrains Mono", monospace;
    font-size: 12px;
    padding: 6px;
    selection-background-color: #1e4620;
}

QLineEdit {
    background-color: #1e293b;
    color: #e2e8f0;
    border: 1px solid #334155;
    border-radius: 4px;
    padding: 6px 10px;
    font-size: 13px;
}

QLineEdit:focus {
    border-color: #4ade80;
}

QPushButton {
    background-color: #1e293b;
    color: #e2e8f0;
    border: 1px solid #334155;
    border-radius: 4px;
    padding: 7px 16px;
    font-weight: 500;
}

QPushButton:hover {
    background-color: #334155;
    border-color: #4ade80;
}

QPushButton:pressed {
    background-color: #0f172a;
}

QPushButton:disabled {
    color: #475569;
    border-color: #1e293b;
}

QPushButton#primary {
    background-color: #16a34a;
    color: white;
    border: none;
}

QPushButton#primary:hover {
    background-color: #15803d;
}

QPushButton#danger {
    background-color: #dc2626;
    color: white;
    border: none;
}

QPushButton#secondary {
    background-color: #2563eb;
    color: white;
    border: none;
}

QPushButton#secondary:hover {
    background-color: #1d4ed8;
}

QComboBox {
    background-color: #1e293b;
    color: #e2e8f0;
    border: 1px solid #334155;
    border-radius: 4px;
    padding: 6px 10px;
}

QComboBox::drop-down {
    border: none;
    padding-right: 8px;
}

QComboBox QAbstractItemView {
    background-color: #1e293b;
    color: #e2e8f0;
    selection-background-color: #16a34a;
}

QSpinBox, QDoubleSpinBox {
    background-color: #1e293b;
    color: #e2e8f0;
    border: 1px solid #334155;
    border-radius: 4px;
    padding: 5px 8px;
}

QTableWidget {
    background-color: #0f172a;
    color: #e2e8f0;
    gridline-color: #1e293b;
    border: 1px solid #1e293b;
    selection-background-color: #1e3a2e;
}

QTableWidget::item {
    padding: 4px 8px;
}

QTableWidget::item:selected {
    background-color: #1e3a2e;
    color: #4ade80;
}

QHeaderView::section {
    background-color: #1e293b;
    color: #94a3b8;
    padding: 6px 8px;
    border: none;
    border-right: 1px solid #0f172a;
    font-weight: 600;
    font-size: 11px;
    text-transform: uppercase;
}

QListWidget {
    background-color: #0f172a;
    color: #e2e8f0;
    border: 1px solid #1e293b;
    border-radius: 4px;
}

QListWidget::item {
    padding: 6px 10px;
}

QListWidget::item:selected {
    background-color: #1e3a2e;
    color: #4ade80;
}

QListWidget::item:hover {
    background-color: #1e293b;
}

QSplitter::handle {
    background-color: #1e293b;
    width: 2px;
    height: 2px;
}

QScrollBar:vertical {
    background-color: #0f172a;
    width: 8px;
    border: none;
}

QScrollBar::handle:vertical {
    background-color: #334155;
    border-radius: 4px;
    min-height: 20px;
}

QScrollBar::handle:vertical:hover {
    background-color: #475569;
}

QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0px;
}

QScrollBar:horizontal {
    background-color: #0f172a;
    height: 8px;
    border: none;
}

QScrollBar::handle:horizontal {
    background-color: #334155;
    border-radius: 4px;
}

QLabel {
    color: #e2e8f0;
}

QLabel#heading {
    color: #4ade80;
    font-size: 16px;
    font-weight: 600;
}

QLabel#subheading {
    color: #94a3b8;
    font-size: 12px;
}

QLabel#code {
    font-family: "JetBrains Mono", monospace;
    color: #4ade80;
    background-color: #020617;
    padding: 4px 8px;
    border-radius: 3px;
}

QGroupBox {
    color: #94a3b8;
    border: 1px solid #1e293b;
    border-radius: 6px;
    margin-top: 10px;
    padding-top: 10px;
    font-weight: 600;
    font-size: 11px;
    text-transform: uppercase;
}

QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 4px;
}

QProgressBar {
    background-color: #1e293b;
    border: none;
    border-radius: 3px;
    height: 6px;
    text-align: center;
}

QProgressBar::chunk {
    background-color: #16a34a;
    border-radius: 3px;
}

QCheckBox {
    spacing: 8px;
}

QCheckBox::indicator {
    width: 16px;
    height: 16px;
    border: 2px solid #334155;
    border-radius: 3px;
    background-color: #1e293b;
}

QCheckBox::indicator:checked {
    background-color: #16a34a;
    border-color: #16a34a;
}

QRadioButton {
    spacing: 8px;
}

QRadioButton::indicator {
    width: 16px;
    height: 16px;
    border: 2px solid #334155;
    border-radius: 8px;
    background-color: #1e293b;
}

QRadioButton::indicator:checked {
    background-color: #16a34a;
    border-color: #16a34a;
}

QStatusBar {
    background-color: #020617;
    color: #4ade80;
    border-top: 1px solid #1e293b;
    font-size: 12px;
}

QMenuBar {
    background-color: #020617;
    color: #e2e8f0;
    border-bottom: 1px solid #1e293b;
}

QMenuBar::item:selected {
    background-color: #1e293b;
}

QMenu {
    background-color: #1e293b;
    color: #e2e8f0;
    border: 1px solid #334155;
    border-radius: 4px;
}

QMenu::item:selected {
    background-color: #16a34a;
    color: white;
}

QToolBar {
    background-color: #0f172a;
    border-bottom: 1px solid #1e293b;
    spacing: 4px;
    padding: 4px;
}

QSplitter {
    background-color: #0f172a;
}

QFrame[frameShape="4"], QFrame[frameShape="5"] {
    background-color: #1e293b;
    border: none;
}
"""
