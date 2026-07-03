"""Vendor submission provider abstraction (issue #10, Phase 1 — credentials only).

Per the #10 design-plan comment, this is the credential/connection-testing
slice: it stores vendor API credentials and can confirm they're present, but
it does **not** submit orders. Both IDT (SciTools Plus API) and Twist
Bioscience (TAPI) gate live API access behind a sales/partnership onboarding
process — there is no public sandbox either vendor exposes without a signed
agreement (checked directly against their own developer pages: IDT requires
booking an integration consultation; Twist requires emailing to register an
IP-whitelisted token). Without a live endpoint to call, `submit()` and
`quote()` deliberately raise rather than fabricate a response — the offline
`synthesis_order.export_order()` file-generation path remains the real
"ordering" feature until a partnership makes live submission possible.

`test_connection()` is honest about this: it reports whether credentials are
saved locally, not a live vendor round-trip that doesn't exist yet.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

from .credentials import has_credentials


class ProviderNotAvailable(RuntimeError):
    """Raised by quote()/submit()/order_status() — no live vendor endpoint yet."""


@dataclass
class ConnectionCheck:
    credentials_present: bool
    live_connection_tested: bool
    message: str


class SubmitProvider(ABC):
    """One implementation per synthesis vendor.

    `vendor_id` matches the key used in `synthesis_order.VENDORS`, so a
    future GUI can offer "test/manage credentials" alongside the existing
    offline export options for the same vendor list.
    """

    vendor_id: str
    required_credential_fields: list[str]
    #: Short, user-facing note on how to obtain this vendor's API access, shown
    #: beside its credential fields in the Settings UI (issue #47).
    credential_help: str = ""
    #: Link to the vendor's own API-access documentation for that note.
    credential_help_url: str = ""

    def test_connection(self) -> ConnectionCheck:
        present = has_credentials(self.vendor_id, self.required_credential_fields)
        if not present:
            return ConnectionCheck(
                credentials_present=False,
                live_connection_tested=False,
                message=(
                    f"No {self.vendor_id} credentials saved yet "
                    f"(needs: {', '.join(self.required_credential_fields)})."
                ),
            )
        return ConnectionCheck(
            credentials_present=True,
            live_connection_tested=False,
            message=(
                f"{self.vendor_id} credentials are saved locally in the OS keychain. "
                "A live connection test isn't available yet — this vendor gates API "
                "access behind a partnership/onboarding process rather than an open "
                "developer sandbox, so there is no endpoint to verify against until "
                "that access is granted."
            ),
        )

    @abstractmethod
    def quote(self, *args, **kwargs):
        raise ProviderNotAvailable(
            f"{self.vendor_id} live ordering is not available — see class docstring "
            "in backend/ordering/providers.py for why. Use the offline export "
            "(backend/formats/synthesis_order.py) instead."
        )

    @abstractmethod
    def submit(self, *args, **kwargs):
        raise ProviderNotAvailable(
            f"{self.vendor_id} live ordering is not available — see class docstring "
            "in backend/ordering/providers.py for why. Use the offline export "
            "(backend/formats/synthesis_order.py) instead."
        )

    @abstractmethod
    def order_status(self, order_id: str):
        raise ProviderNotAvailable(
            f"{self.vendor_id} live ordering is not available — see class docstring "
            "in backend/ordering/providers.py for why."
        )


class IDTProvider(SubmitProvider):
    """IDT SciTools Plus API — access requires booking an integration
    consultation with IDT (idtdna.com/pages/products/gmp-oem-and-integrations/
    integrations/scitools-plus-api); exact credential shape (API key vs.
    OAuth client id/secret) should be confirmed once that access exists.
    """

    vendor_id = "idt"
    required_credential_fields = ["api_key"]
    credential_help = (
        "IDT's SciTools Plus API isn't self-serve — access is granted by booking "
        "an integration consultation with IDT, who then provide your API "
        "credentials. There's no public sandbox key to generate yourself."
    )
    credential_help_url = (
        "https://www.idtdna.com/pages/products/gmp-oem-and-integrations/"
        "integrations/scitools-plus-api"
    )

    def quote(self, *args, **kwargs):
        return super().quote(*args, **kwargs)

    def submit(self, *args, **kwargs):
        return super().submit(*args, **kwargs)

    def order_status(self, order_id: str):
        return super().order_status(order_id)


class TwistProvider(SubmitProvider):
    """Twist Bioscience TAPI — access requires emailing Twist to register an
    account email + whitelisted IP address, after which a one-time-use
    token link is issued (twistbioscience.com/tapi).
    """

    vendor_id = "twist"
    required_credential_fields = ["api_token", "account_email"]
    credential_help = (
        "Twist's TAPI access is arranged by emailing Twist to register your "
        "account email and a whitelisted IP address; they then issue a "
        "one-time-use token link. Use that account email and token here."
    )
    credential_help_url = "https://www.twistbioscience.com/tapi"

    def quote(self, *args, **kwargs):
        return super().quote(*args, **kwargs)

    def submit(self, *args, **kwargs):
        return super().submit(*args, **kwargs)

    def order_status(self, order_id: str):
        return super().order_status(order_id)


PROVIDERS: dict[str, type[SubmitProvider]] = {
    "idt": IDTProvider,
    "twist": TwistProvider,
}
