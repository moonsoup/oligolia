"""#61: ClinVar significance was taken from a gene search and applied to everything.

Three defects, all in `backend/routers/variants.py`:

  1. the query was `{gene}[Gene Name]` and the answer was `records[0]`, so every
     variant in BRCA1 got whatever the first returned BRCA1 record said —
     regardless of position or allele;
  2. the enum was matched by SUBSTRING, so "Likely pathogenic" matched PATHOGENIC
     (declared first, and "pathogenic" is a substring of it), and ClinVar's
     current "Conflicting classifications of pathogenicity" did too;
  3. `except Exception: pass` hid every failure.

The matcher is pure, so all of this is testable without touching the network.
Fixtures use the real esummary shape: `variation_set[].variation_loc[]` with
assembly/chr/start/ref/alt, plus `germline_classification` (current) and
`clinical_significance` (legacy).
"""

from __future__ import annotations

import pytest

from backend.models.variant import ClinicalSignificance, Variant, VariantType
from backend.routers.variants import (
    classify_significance,
    matching_clinvar_records,
    significance_for_variant,
)


def _record(chrom: str, pos: int, ref: str, alt: str, description: str,
            *, assembly: str = "GRCh38", legacy: bool = False) -> dict:
    key = "clinical_significance" if legacy else "germline_classification"
    return {
        "obj_type": "single nucleotide variant",
        key: {"description": description, "last_evaluated": "2024/01/01", "review_status": "criteria provided"},
        "variation_set": [{
            "variation_name": f"NM_007294.4(BRCA1):c.{pos}{ref}>{alt}",
            "canonical_spdi": f"NC_000017.11:{pos - 1}:{ref}:{alt}",
            "variation_loc": [{
                "assembly_name": assembly, "chr": chrom,
                "start": str(pos), "stop": str(pos),
                "ref": ref, "alt": alt, "status": "current",
            }],
        }],
    }


def _variant(chrom="17", pos=43045703, ref="A", alt=("G",), gene="BRCA1") -> Variant:
    return Variant(chrom=chrom, pos=pos, ref=ref, alt=list(alt),
                   variant_type=VariantType.SNP, gene=gene)


# ── classify_significance: exact, not substring ──────────────────────────────

@pytest.mark.parametrize("description,expected", [
    ("Pathogenic", ClinicalSignificance.PATHOGENIC),
    ("Likely pathogenic", ClinicalSignificance.LIKELY_PATHOGENIC),
    ("Benign", ClinicalSignificance.BENIGN),
    ("Likely benign", ClinicalSignificance.LIKELY_BENIGN),
    ("Uncertain significance", ClinicalSignificance.VUS),
    ("Conflicting classifications of pathogenicity", ClinicalSignificance.CONFLICTING),
    ("Not provided", ClinicalSignificance.NOT_PROVIDED),
])
def test_known_classifications_map_exactly(description: str, expected: ClinicalSignificance) -> None:
    assert classify_significance(description) is expected, description


def test_likely_pathogenic_is_not_pathogenic() -> None:
    """The substring bug: 'pathogenic' is inside 'likely pathogenic'."""
    assert classify_significance("Likely pathogenic") is ClinicalSignificance.LIKELY_PATHOGENIC
    assert classify_significance("Likely pathogenic") is not ClinicalSignificance.PATHOGENIC


def test_conflicting_classifications_is_not_pathogenic() -> None:
    """ClinVar's current wording contains 'pathogenicity'. It is not a verdict."""
    got = classify_significance("Conflicting classifications of pathogenicity")
    assert got is not ClinicalSignificance.PATHOGENIC, got
    # The enum already had the right member for this, which the substring match
    # could never reach because PATHOGENIC is declared first.
    assert got is ClinicalSignificance.CONFLICTING, got


def test_an_unrecognised_description_is_not_guessed() -> None:
    assert classify_significance("something ClinVar invented last week") is None
    assert classify_significance("") is None
    # "not provided" IS a real ClinVar value with its own enum member.
    assert classify_significance("not provided") is ClinicalSignificance.NOT_PROVIDED


def test_classification_is_case_and_space_insensitive() -> None:
    assert classify_significance("  PATHOGENIC  ") is ClinicalSignificance.PATHOGENIC
    assert classify_significance("likely  benign") is ClinicalSignificance.LIKELY_BENIGN


# ── matching_clinvar_records: position and allele must agree ─────────────────

def test_a_record_at_the_same_position_and_allele_matches() -> None:
    v = _variant()
    recs = [_record("17", 43045703, "A", "G", "Pathogenic")]
    assert matching_clinvar_records(v, recs) == recs


def test_a_record_at_a_different_position_does_not_match() -> None:
    """The heart of #61: a neighbour's verdict is not this variant's verdict."""
    v = _variant(pos=43045703)
    recs = [_record("17", 43099999, "A", "G", "Pathogenic")]
    assert matching_clinvar_records(v, recs) == []


def test_a_record_with_a_different_alt_does_not_match() -> None:
    v = _variant(alt=("G",))
    recs = [_record("17", 43045703, "A", "T", "Pathogenic")]
    assert matching_clinvar_records(v, recs) == []


def test_a_record_on_a_different_chromosome_does_not_match() -> None:
    v = _variant(chrom="17")
    recs = [_record("13", 43045703, "A", "G", "Pathogenic")]
    assert matching_clinvar_records(v, recs) == []


def test_the_chr_prefix_is_tolerated_on_either_side() -> None:
    """VCFs write both '17' and 'chr17'; ClinVar writes '17'."""
    recs = [_record("17", 43045703, "A", "G", "Pathogenic")]
    assert matching_clinvar_records(_variant(chrom="chr17"), recs) == recs
    assert matching_clinvar_records(_variant(chrom="17"), recs) == recs


def test_any_of_a_multi_alt_variant_may_match() -> None:
    v = _variant(alt=("T", "G"))
    recs = [_record("17", 43045703, "A", "G", "Pathogenic")]
    assert matching_clinvar_records(v, recs) == recs


def test_only_the_matching_record_is_kept_from_a_gene_wide_result() -> None:
    """Exactly the #61 scenario: five BRCA1 records, one of them ours."""
    v = _variant(pos=43045703, alt=("G",))
    recs = [
        _record("17", 43044295, "T", "C", "Pathogenic"),
        _record("17", 43045703, "A", "G", "Benign"),
        _record("17", 43047643, "G", "A", "Pathogenic"),
    ]
    matched = matching_clinvar_records(v, recs)
    assert len(matched) == 1
    assert matched[0] is recs[1]


# ── significance_for_variant: the two together ───────────────────────────────

def test_the_variants_own_significance_is_returned() -> None:
    v = _variant(pos=43045703, alt=("G",))
    recs = [
        _record("17", 43044295, "T", "C", "Pathogenic"),
        _record("17", 43045703, "A", "G", "Benign"),
    ]
    assert significance_for_variant(v, recs) is ClinicalSignificance.BENIGN


def test_no_matching_record_returns_nothing_not_a_neighbour() -> None:
    """#61 in one assertion."""
    v = _variant(pos=43045703)
    recs = [_record("17", 43044295, "T", "C", "Pathogenic")]
    assert significance_for_variant(v, recs) is None


def test_an_empty_result_returns_nothing() -> None:
    assert significance_for_variant(_variant(), []) is None


def test_the_legacy_clinical_significance_field_is_still_read() -> None:
    """Older esummary payloads use clinical_significance, not germline_classification."""
    v = _variant(pos=43045703, alt=("G",))
    recs = [_record("17", 43045703, "A", "G", "Likely benign", legacy=True)]
    assert significance_for_variant(v, recs) is ClinicalSignificance.LIKELY_BENIGN


def test_a_conflicting_record_is_reported_as_conflicting_not_pathogenic() -> None:
    v = _variant(pos=43045703, alt=("G",))
    recs = [_record("17", 43045703, "A", "G", "Conflicting classifications of pathogenicity")]
    got = significance_for_variant(v, recs)
    assert got is not ClinicalSignificance.PATHOGENIC, got
    assert got is ClinicalSignificance.CONFLICTING, got


def test_a_non_grch38_location_is_not_matched_against_grch38_coordinates() -> None:
    """GRCh37 and GRCh38 coordinates differ by megabases in places."""
    v = _variant(pos=43045703)
    recs = [_record("17", 43045703, "A", "G", "Pathogenic", assembly="GRCh37")]
    assert matching_clinvar_records(v, recs) == []
