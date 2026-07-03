"""Settings panel — vendor ordering credentials (issue #10, Phase 1 GUI).

A UI over the existing backend/ordering credential plumbing: enter, test, and
clear per-vendor API credentials. Credentials go to the OS keychain via
`keyring` (never plaintext), and NOTHING is sent to any vendor — live
submission isn't available yet (both IDT and Twist gate API access behind a
partnership). "Test connection" reports whether credentials are saved
locally, not a live round-trip.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QLineEdit,
    QPushButton, QGroupBox, QMessageBox, QScrollArea,
)

from backend.formats.synthesis_order import VENDORS
from backend.ordering import (
    PROVIDERS, delete_credential, get_credential, has_credentials, save_credential,
)

_SECRET_HINTS = ("key", "token", "secret", "password")


def _is_secret(field: str) -> bool:
    return any(h in field.lower() for h in _SECRET_HINTS)


def _vendor_label(vendor_id: str) -> str:
    profile = VENDORS.get(vendor_id)
    return profile.display_name if profile else vendor_id.upper()


class _VendorCredentials(QGroupBox):
    """One vendor's credential fields + Save / Test / Clear."""

    def __init__(self, vendor_id: str) -> None:
        super().__init__(_vendor_label(vendor_id))
        self._vendor_id = vendor_id
        self._provider = PROVIDERS[vendor_id]()
        self._fields: dict[str, QLineEdit] = {}

        layout = QVBoxLayout(self)
        form = QFormLayout()
        for field in self._provider.required_credential_fields:
            edit = QLineEdit()
            if _is_secret(field):
                edit.setEchoMode(QLineEdit.EchoMode.Password)
            self._fields[field] = edit
            form.addRow(field.replace("_", " ").title() + ":", edit)
        layout.addLayout(form)

        # How to obtain this vendor's key (issue #47) — a short note + a link to
        # the vendor's own API-access docs, so a user isn't left guessing where
        # a key comes from (neither vendor offers a self-serve key).
        help_text = getattr(self._provider, "credential_help", "")
        help_url = getattr(self._provider, "credential_help_url", "")
        if help_text or help_url:
            parts = []
            if help_text:
                parts.append(help_text)
            if help_url:
                parts.append(
                    f'<a href="{help_url}">How to get {_vendor_label(vendor_id)} '
                    f"API access &rarr;</a>"
                )
            help_label = QLabel(" ".join(parts))
            help_label.setWordWrap(True)
            help_label.setObjectName("subheading")
            help_label.setOpenExternalLinks(True)
            help_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextBrowserInteraction)
            layout.addWidget(help_label)

        btn_row = QHBoxLayout()
        btn_save = QPushButton("Save credentials")
        btn_save.setObjectName("primary")
        btn_save.clicked.connect(self._save)
        btn_test = QPushButton("Test connection")
        btn_test.clicked.connect(self._test)
        btn_clear = QPushButton("Clear")
        btn_clear.setObjectName("danger")
        btn_clear.clicked.connect(self._clear)
        btn_row.addWidget(btn_save)
        btn_row.addWidget(btn_test)
        btn_row.addWidget(btn_clear)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        self._status = QLabel("")
        self._status.setObjectName("subheading")
        self._status.setWordWrap(True)
        layout.addWidget(self._status)

        self._refresh()

    def _refresh(self) -> None:
        """Reflect stored-credential state without ever showing secret values."""
        try:
            for field, edit in self._fields.items():
                stored = get_credential(self._vendor_id, field)
                if _is_secret(field):
                    edit.clear()
                    edit.setPlaceholderText("•••••• saved — leave blank to keep" if stored
                                            else "not set")
                else:
                    # Non-secret fields (e.g. account email) can be shown/edited.
                    edit.setText(stored or "")
                    edit.setPlaceholderText("not set")
            present = has_credentials(self._vendor_id, self._provider.required_credential_fields)
            self._status.setText("✓ credentials saved (OS keychain)" if present
                                 else "No credentials saved yet.")
        except Exception as e:  # no OS keychain backend available, etc.
            self._status.setText(f"Keychain unavailable: {e}")

    def _save(self) -> None:
        try:
            saved_any = False
            for field, edit in self._fields.items():
                val = edit.text().strip()
                if val:  # blank => keep any existing value (esp. for secrets)
                    save_credential(self._vendor_id, field, val)
                    saved_any = True
            self._refresh()
            self._status.setText("Saved to OS keychain." if saved_any
                                 else "Nothing to save (all fields blank).")
        except Exception as e:
            QMessageBox.critical(self, "Save failed", f"Could not store credentials:\n{e}")

    def _test(self) -> None:
        # No network — reports whether credentials are present locally.
        try:
            check = self._provider.test_connection()
        except Exception as e:
            self._status.setText(f"Test failed: {e}")
            return
        mark = "✓" if check.credentials_present else "•"
        self._status.setText(f"{mark} {check.message}")

    def _clear(self) -> None:
        if QMessageBox.question(
            self, "Clear credentials",
            f"Remove saved {_vendor_label(self._vendor_id)} credentials from the keychain?",
        ) != QMessageBox.StandardButton.Yes:
            return
        try:
            for field in self._provider.required_credential_fields:
                delete_credential(self._vendor_id, field)
        except Exception as e:
            QMessageBox.critical(self, "Clear failed", str(e))
        self._refresh()


class SettingsPanel(QWidget):
    def __init__(self) -> None:
        super().__init__()
        outer = QVBoxLayout(self)
        outer.setContentsMargins(8, 8, 8, 8)

        heading = QLabel("Vendor Ordering Credentials")
        heading.setObjectName("heading")
        outer.addWidget(heading)

        intro = QLabel(
            "Store API credentials for synthesis vendors. Credentials are kept in "
            "your operating system's keychain (never in a plaintext file), and "
            "<b>nothing is sent to any vendor from here</b> — live ordering isn't "
            "available yet, so \"Test connection\" only confirms your credentials "
            "are saved locally. The current way to order is still the offline "
            "<i>Export for Synthesis</i> file you upload on the vendor's own portal."
        )
        intro.setWordWrap(True)
        intro.setObjectName("subheading")
        outer.addWidget(intro)

        # Scrollable vendor list (room for more vendors/settings later).
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        container = QWidget()
        col = QVBoxLayout(container)
        for vendor_id in PROVIDERS:
            col.addWidget(_VendorCredentials(vendor_id))
        col.addStretch()
        scroll.setWidget(container)
        outer.addWidget(scroll)
