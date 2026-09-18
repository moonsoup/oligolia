"""#67.1: symbolic ALT alleles were typed by string length.

`_infer_type("A", ["<DEL>"])` compared len("A") with len("<DEL>") and returned
INS — a structural deletion reported as an insertion. Every symbolic form
(`<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, breakends) went through the same length
heuristic, so all of them were mislabelled, silently.

VCF 4.x defines these explicitly, so they are matched explicitly now, and an
unrecognised symbolic allele is SV rather than a guess from its spelling.
"""

from __future__ import annotations

import pytest

from backend.formats.vcf import _infer_type
from backend.models.variant import VariantType


@pytest.mark.parametrize("alt,expected", [
    ("<DEL>", VariantType.DEL),
    ("<INS>", VariantType.INS),
    ("<DUP>", VariantType.CNV),
    ("<CNV>", VariantType.CNV),
    ("<INV>", VariantType.SV),
    ("<DUP:TANDEM>", VariantType.CNV),
    ("<DEL:ME:ALU>", VariantType.DEL),
    ("<INS:ME:L1>", VariantType.INS),
])
def test_symbolic_alts_are_typed_by_their_symbol(alt: str, expected: VariantType) -> None:
    assert _infer_type("A", [alt]) is expected, alt


def test_an_unknown_symbolic_alt_is_sv_not_a_length_guess() -> None:
    """The honest answer for something VCF does not define: structural, unspecified."""
    assert _infer_type("A", ["<SOMETHING_NEW>"]) is VariantType.SV
    assert _infer_type("ACGTACGT", ["<X>"]) is VariantType.SV


@pytest.mark.parametrize("alt", ["]17:198982]A", "A]17:198982]", "[13:123457[A", "A[13:123457["])
def test_breakend_notation_is_sv(alt: str) -> None:
    """VCF 4.x breakends. Length comparison made these arbitrary."""
    assert _infer_type("A", [alt]) is VariantType.SV, alt


def test_the_plain_cases_are_unchanged() -> None:
    """Regression guard: the length heuristic was right for literal alleles."""
    assert _infer_type("A", ["G"]) is VariantType.SNP
    assert _infer_type("AC", ["GT"]) is VariantType.SNP          # MNP
    assert _infer_type("A", ["ACGT"]) is VariantType.INS
    assert _infer_type("ACGT", ["A"]) is VariantType.DEL
    assert _infer_type("A", ["."]) is VariantType.UNKNOWN
    assert _infer_type("A", ["*"]) is VariantType.UNKNOWN


def test_a_symbolic_alt_beside_a_literal_one_is_still_recognised() -> None:
    """Multi-ALT records are real; the symbolic one must not be skipped over."""
    assert _infer_type("A", ["<DEL>", "G"]) is VariantType.DEL
    assert _infer_type("A", [".", "<DUP>"]) is VariantType.CNV


def test_symbolic_alts_parse_end_to_end() -> None:
    """Through the real parser, not just the helper."""
    from backend.formats.vcf import parse_vcf

    vcf = (
        "##fileformat=VCFv4.2\n"
        "##ALT=<ID=DEL,Description=\"Deletion\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t1000\tsv1\tA\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=2000\n"
        "chr1\t3000\tsv2\tA\t<DUP>\t50\tPASS\tSVTYPE=DUP;END=4000\n"
        "chr1\t5000\trs1\tA\tG\t50\tPASS\t.\n"
    )
    variants = parse_vcf(vcf)
    assert [v.variant_type for v in variants] == [
        VariantType.DEL, VariantType.CNV, VariantType.SNP,
    ], [v.variant_type for v in variants]
