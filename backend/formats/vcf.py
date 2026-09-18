"""VCF (Variant Call Format) parser — pure Python, no pysam dependency required."""

from __future__ import annotations
from typing import Iterator, TextIO
from ..models.variant import Variant, VariantType


#: Symbolic ALT alleles from the VCF 4.x spec, mapped by their leading keyword so
#: that subtypes like <DEL:ME:ALU> and <DUP:TANDEM> resolve too.
#:
#: These must be matched BEFORE any length comparison. `len("<DEL>") > len("A")`,
#: so the old length heuristic reported every structural deletion as an insertion
#: (#67). A symbolic allele's spelling says nothing about the variant's size.
_SYMBOLIC_ALT = {
    "DEL": VariantType.DEL,
    "INS": VariantType.INS,
    "DUP": VariantType.CNV,   # a duplication is a copy-number gain
    "CNV": VariantType.CNV,
    "INV": VariantType.SV,
    "TRA": VariantType.SV,
    "BND": VariantType.SV,
}


def _classify_alt(ref: str, alt: str) -> VariantType | None:
    """Type one ALT allele, or None if it carries no information."""
    if alt in (".", "*", ""):
        return None

    # Breakend notation: A[13:123457[ , ]17:198982]A and the other two forms.
    if "[" in alt or "]" in alt:
        return VariantType.SV

    # Symbolic allele: <DEL>, <DUP:TANDEM>, <DEL:ME:ALU>, ...
    if alt.startswith("<") and alt.endswith(">"):
        keyword = alt[1:-1].split(":", 1)[0].upper()
        return _SYMBOLIC_ALT.get(keyword, VariantType.SV)

    # Literal alleles: the length comparison, which was always right for these.
    if len(ref) == len(alt):
        return VariantType.SNP  # SNP, or an MNP treated as the SNP class
    return VariantType.INS if len(alt) > len(ref) else VariantType.DEL


def _infer_type(ref: str, alts: list[str]) -> VariantType:
    for alt in alts:
        kind = _classify_alt(ref, alt)
        if kind is not None:
            return kind
    return VariantType.UNKNOWN


def _parse_info(info_str: str) -> dict:
    if info_str == ".":
        return {}
    result: dict = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            result[k] = v
        else:
            result[item] = True
    return result


def parse_vcf(source: str | TextIO) -> list[Variant]:
    if isinstance(source, str):
        lines = source.splitlines()
    else:
        lines = list(source)

    variants = []
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        chrom, pos_str, id_, ref, alt_str = parts[:5]
        qual_str = parts[5] if len(parts) > 5 else "."
        filter_str = parts[6] if len(parts) > 6 else "."
        info_str = parts[7] if len(parts) > 7 else "."

        alts = [a for a in alt_str.split(",") if a not in (".", "")]
        qual = None if qual_str in (".", "") else float(qual_str)
        filters = [] if filter_str in (".", "PASS", "") else filter_str.split(";")
        info = _parse_info(info_str)

        variants.append(Variant(
            chrom=chrom,
            pos=int(pos_str),
            ref=ref,
            alt=alts,
            id=id_,
            qual=qual,
            filter=filters,
            info=info,
            variant_type=_infer_type(ref, alts),
            gene=info.get("GENEINFO", "").split(":")[0] if "GENEINFO" in info else "",
        ))
    return variants


def write_vcf(variants: list[Variant], source_name: str = "oligolia") -> str:
    lines = [
        "##fileformat=VCFv4.3",
        f"##source={source_name}",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for v in variants:
        alt = ",".join(v.alt) if v.alt else "."
        qual = str(v.qual) if v.qual is not None else "."
        filter_ = ";".join(v.filter) if v.filter else "PASS"
        info_parts = [f"{k}={val}" if val is not True else k
                      for k, val in v.info.items()]
        info = ";".join(info_parts) if info_parts else "."
        lines.append(f"{v.chrom}\t{v.pos}\t{v.id}\t{v.ref}\t{alt}\t{qual}\t{filter_}\t{info}")
    return "\n".join(lines) + "\n"


def parse_vcf_iter(source: TextIO) -> Iterator[Variant]:
    """Memory-efficient streaming parser for large VCF files."""
    for line in source:
        line = line.rstrip("\n")
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        chrom, pos_str, id_, ref, alt_str = parts[:5]
        alts = [a for a in alt_str.split(",") if a not in (".", "")]
        yield Variant(
            chrom=chrom, pos=int(pos_str), ref=ref, alt=alts, id=id_,
            variant_type=_infer_type(ref, alts),
        )
