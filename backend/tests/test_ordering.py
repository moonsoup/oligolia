"""Tests for backend/ordering/ — credential storage and provider abstraction.

Uses an in-memory fake keyring backend so these tests never touch the real
OS keychain (which would hang on a permission prompt in CI/automation) and
never leave real secrets behind on the test machine.
"""
from __future__ import annotations

import keyring
import pytest
from keyring.backend import KeyringBackend
from keyring.errors import PasswordDeleteError

from backend.ordering import credentials
from backend.ordering.providers import (
    IDTProvider, TwistProvider, ProviderNotAvailable,
)


class _InMemoryKeyring(KeyringBackend):
    priority = 1  # highest priority so keyring picks this backend

    def __init__(self) -> None:
        self._store: dict[tuple[str, str], str] = {}

    def get_password(self, service, username):
        return self._store.get((service, username))

    def set_password(self, service, username, password):
        self._store[(service, username)] = password

    def delete_password(self, service, username):
        if (service, username) not in self._store:
            raise PasswordDeleteError("not found")
        del self._store[(service, username)]


@pytest.fixture(autouse=True)
def fake_keyring():
    original = keyring.get_keyring()
    fake = _InMemoryKeyring()
    keyring.set_keyring(fake)
    yield fake
    keyring.set_keyring(original)


# ── credentials.py ──────────────────────────────────────────────────────────

def test_save_and_get_credential():
    credentials.save_credential("idt", "api_key", "secret-123")
    assert credentials.get_credential("idt", "api_key") == "secret-123"


def test_get_missing_credential_returns_none():
    assert credentials.get_credential("idt", "api_key") is None


def test_delete_credential():
    credentials.save_credential("twist", "api_token", "tok-abc")
    assert credentials.delete_credential("twist", "api_token") is True
    assert credentials.get_credential("twist", "api_token") is None


def test_delete_missing_credential_returns_false():
    assert credentials.delete_credential("twist", "nope") is False


def test_has_credentials_all_present():
    credentials.save_credential("twist", "api_token", "tok")
    credentials.save_credential("twist", "account_email", "a@b.com")
    assert credentials.has_credentials("twist", ["api_token", "account_email"]) is True


def test_has_credentials_missing_one():
    credentials.save_credential("twist", "api_token", "tok")
    assert credentials.has_credentials("twist", ["api_token", "account_email"]) is False


def test_credentials_are_namespaced_per_vendor():
    credentials.save_credential("idt", "api_key", "idt-key")
    credentials.save_credential("twist", "api_key", "twist-key")
    assert credentials.get_credential("idt", "api_key") == "idt-key"
    assert credentials.get_credential("twist", "api_key") == "twist-key"


# ── providers.py ─────────────────────────────────────────────────────────────

def test_idt_test_connection_no_credentials():
    result = IDTProvider().test_connection()
    assert result.credentials_present is False
    assert result.live_connection_tested is False
    assert "no idt credentials" in result.message.lower()


def test_idt_test_connection_with_credentials():
    credentials.save_credential("idt", "api_key", "k")
    result = IDTProvider().test_connection()
    assert result.credentials_present is True
    assert result.live_connection_tested is False
    assert "saved locally" in result.message.lower()


def test_twist_requires_both_fields():
    credentials.save_credential("twist", "api_token", "tok")
    # account_email still missing
    result = TwistProvider().test_connection()
    assert result.credentials_present is False


def test_twist_test_connection_with_both_fields():
    credentials.save_credential("twist", "api_token", "tok")
    credentials.save_credential("twist", "account_email", "a@b.com")
    result = TwistProvider().test_connection()
    assert result.credentials_present is True


def test_submit_raises_provider_not_available():
    with pytest.raises(ProviderNotAvailable):
        IDTProvider().submit()


def test_quote_raises_provider_not_available():
    with pytest.raises(ProviderNotAvailable):
        TwistProvider().quote()


def test_order_status_raises_provider_not_available():
    with pytest.raises(ProviderNotAvailable):
        IDTProvider().order_status("some-id")
