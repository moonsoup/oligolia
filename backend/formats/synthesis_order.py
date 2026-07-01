"""
Synthesis-order export — builds the upload file a synthesis vendor's own
ordering portal expects for a sequence (or sub-range).

This only generates a file the user uploads themselves on the vendor's site;
it never submits an order. Column layouts are taken from each vendor's
published ordering guide where noted; vendors marked verified=False have no
confirmed public spec and should be double-checked against the vendor's own
template before use.
"""

from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO, StringIO

from openpyxl import Workbook

from ..models.sequence import Sequence


@dataclass(frozen=True)
class VendorProfile:
    key: str
    display_name: str
    file_format: str  # "xlsx" | "csv"
    portal_url: str
    requires_account: bool
    verified: bool  # export format confirmed against the vendor's own published spec
    instructions: list[str]


VENDORS: dict[str, VendorProfile] = {
    "genewiz": VendorProfile(
        key="genewiz",
        display_name="GENEWIZ from Azenta (Gene Synthesis)",
        file_format="xlsx",
        portal_url="https://www.genewiz.com",
        requires_account=True,
        verified=True,
        instructions=[
            "Log into your GENEWIZ/Azenta account.",
            "Select the ‘Gene Synthesis’ tab and choose a service (e.g. PriorityGENE).",
            "Switch the order form to ‘Grid View’ (top right corner).",
            "Click Download/Upload → Upload Excel, and select the file generated here.",
            "Review vector, cloning, and DNA-prep options, then submit for a quote.",
        ],
    ),
    "idt": VendorProfile(
        key="idt",
        display_name="Integrated DNA Technologies (gBlocks/eBlocks)",
        file_format="csv",
        portal_url="https://www.idtdna.com/site/order/gblockentry",
        requires_account=False,
        verified=True,
        instructions=[
            "Open IDT’s gBlocks Gene Fragments order entry tool (no login needed to start).",
            "Use the bulk-paste option and choose ‘Comma-separated’ as the delimiter.",
            "Paste the contents of the generated CSV (Name,Sequence).",
            "Review IDT’s complexity/length warnings, then add to cart.",
        ],
    ),
    "twist": VendorProfile(
        key="twist",
        display_name="Twist Bioscience (Gene Fragments)",
        file_format="csv",
        portal_url="https://www.twistbioscience.com/twist-ordering-platform",
        requires_account=True,
        verified=True,
        instructions=[
            "Log into the Twist Ordering Platform.",
            "Start a new Gene Fragment order and use the table-import option.",
            "Upload the generated CSV (Name,Sequence) — Twist also accepts XLSX/TSV/FASTA.",
            "Review synthesis-feasibility flags, then submit for pricing.",
        ],
    ),
    "eurofins": VendorProfile(
        key="eurofins",
        display_name="Eurofins Genomics (unverified format)",
        file_format="csv",
        portal_url="https://eurofinsgenomics.com/en/orderpages/oligos/custom-dna-oligos-in-plates-order-page/",
        requires_account=True,
        verified=False,
        instructions=[
            "Eurofins uses its own plate-order Excel template (named tabs like "
            "‘Sequence Info A1–A12’); its exact column layout wasn't confirmed from public docs.",
            "This export is a generic Name,Sequence CSV — cross-check it against "
            "Eurofins’ current template before uploading.",
            "Log into your Eurofins Genomics account and download their plate template to compare columns.",
        ],
    ),
}


def export_order(seq: Sequence, vendor: str, start: int = 0, end: int | None = None) -> tuple[bytes, str, list[str]]:
    """Build a vendor-formatted order file for a sequence or sub-range.

    Returns (file_bytes, suggested_filename, hand-off instructions). Does not
    submit anything — the caller writes the bytes to disk for the user to
    upload on the vendor's own portal.
    """
    profile = VENDORS.get(vendor)
    if profile is None:
        raise ValueError(f"Unknown vendor: {vendor}")

    end = len(seq.seq) if end is None else end
    sub_seq = seq.seq[start:end]
    if not sub_seq:
        raise ValueError("Selected range is empty")
    name = seq.name or seq.id

    if profile.file_format == "xlsx":
        data = _genewiz_xlsx(name, sub_seq)
    else:
        data = _name_sequence_csv(name, sub_seq)

    filename = f"{name}_{vendor}_order.{profile.file_format}"
    return data, filename, profile.instructions


def _genewiz_xlsx(name: str, sequence: str) -> bytes:
    """GENEWIZ Grid View columns, per their Gene Synthesis Ordering Guide."""
    wb = Workbook()
    ws = wb.active
    ws.title = "Sequences"
    ws.append(["Sequence Name", "Sequence Type", "5' Flanking", "Sequence", "3' Flanking"])
    ws.append([name, "DNA", "", sequence, ""])
    buf = BytesIO()
    wb.save(buf)
    return buf.getvalue()


def _name_sequence_csv(name: str, sequence: str) -> bytes:
    buf = StringIO()
    buf.write("Name,Sequence\n")
    buf.write(f"{name},{sequence}\n")
    return buf.getvalue().encode("utf-8")
