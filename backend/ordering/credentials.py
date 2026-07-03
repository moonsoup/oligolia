"""Vendor API credential storage — OS keychain via `keyring`, never plaintext.

Per the #10 ordering-integration design plan: QSettings is plaintext on disk
and unsuitable for API tokens. This module is the only place vendor
credentials touch disk, and it never writes them itself — `keyring` hands
storage off to the platform's real secret store (macOS Keychain, Windows
Credential Locker, or a Secret Service provider on Linux).

This module has zero network access and works with zero vendor credentials
present — that's the point: it can be built, tested, and shipped before any
vendor partnership exists, de-risking the auth/storage piece independently
of whether live submission ever ships.
"""
from __future__ import annotations

import os

import keyring
from keyring.errors import PasswordDeleteError

_DEFAULT_SERVICE_PREFIX = "oligolia-ordering"


def _service_prefix() -> str:
    """Keychain namespace prefix, overridable via ``OLIGOLIA_KEYCHAIN_NAMESPACE``.

    A test harness that drives this app in-process (e.g. ProjectMan's playtest
    tool) shares this process's real OS-keychain access with no isolation, so it
    can point at an isolated, always-empty namespace instead of touching a real
    user's saved credentials. Read per call (not cached) so setting the env var
    takes effect for subsequent calls; defaults to the real value, so a normal
    installed copy of the app is unaffected.
    """
    return os.environ.get("OLIGOLIA_KEYCHAIN_NAMESPACE", _DEFAULT_SERVICE_PREFIX)


def _service_name(vendor: str) -> str:
    return f"{_service_prefix()}:{vendor.lower()}"


def save_credential(vendor: str, field: str, value: str) -> None:
    """Store one credential field (e.g. 'api_key') for a vendor."""
    keyring.set_password(_service_name(vendor), field, value)


def get_credential(vendor: str, field: str) -> str | None:
    """Return a stored credential field, or None if never set."""
    return keyring.get_password(_service_name(vendor), field)


def delete_credential(vendor: str, field: str) -> bool:
    """Remove one stored credential field. Returns False if none existed."""
    try:
        keyring.delete_password(_service_name(vendor), field)
        return True
    except PasswordDeleteError:
        return False


def has_credentials(vendor: str, fields: list[str]) -> bool:
    """True only if every required field is present for this vendor."""
    return all(get_credential(vendor, f) is not None for f in fields)
