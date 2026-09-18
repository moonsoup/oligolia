"""Variant annotation and analysis endpoints."""

from fastapi import APIRouter, HTTPException
from ..models.variant import Variant, VariantAnnotationRequest, VariantAnnotationResponse, ClinicalSignificance
from ..services import NCBIClient, GnomADClient

router = APIRouter(prefix="/variants", tags=["variants"])

_ncbi = NCBIClient()
_gnomad = GnomADClient()

#: ClinVar classification strings mapped EXACTLY, never by substring.
#:
#: The previous code iterated the enum and took the first whose value appeared
#: anywhere in the description. `ClinicalSignificance.PATHOGENIC` is declared
#: first and its value "Pathogenic" is a substring of "Likely pathogenic", so
#: every likely-pathogenic variant was reported as pathogenic — and so was
#: ClinVar's current "Conflicting classifications of pathogenicity", which is not
#: a verdict at all (#61).
#:
#: Keys are normalised: lowercased, with runs of whitespace collapsed.
CLINVAR_CLASSIFICATIONS: dict[str, ClinicalSignificance] = {
    "pathogenic": ClinicalSignificance.PATHOGENIC,
    "likely pathogenic": ClinicalSignificance.LIKELY_PATHOGENIC,
    "pathogenic/likely pathogenic": ClinicalSignificance.LIKELY_PATHOGENIC,
    "benign": ClinicalSignificance.BENIGN,
    "likely benign": ClinicalSignificance.LIKELY_BENIGN,
    "benign/likely benign": ClinicalSignificance.LIKELY_BENIGN,
    "uncertain significance": ClinicalSignificance.VUS,
    "uncertain risk allele": ClinicalSignificance.VUS,
    # Explicitly NOT a verdict, and the enum already had the right member for it.
    # The substring match reported these as PATHOGENIC.
    "conflicting classifications of pathogenicity": ClinicalSignificance.CONFLICTING,
    "conflicting interpretations of pathogenicity": ClinicalSignificance.CONFLICTING,
    "not provided": ClinicalSignificance.NOT_PROVIDED,
    "no classification provided": ClinicalSignificance.NOT_PROVIDED,
}

#: The assembly our variant coordinates are assumed to be on. GRCh37 and GRCh38
#: positions differ by megabases in places, so a location on the other assembly
#: is not a match.
CLINVAR_ASSEMBLY = "GRCh38"

#: How many gene records to pull before filtering to the variant. The old code
#: asked for 5 and then took the first, so for any gene with more than a handful
#: of submissions the right record was usually not even fetched.
CLINVAR_MAX_RECORDS = 200


def classify_significance(description: str) -> ClinicalSignificance | None:
    """Map a ClinVar classification string, or None if it is not a known one.

    Returning None is the point: an unrecognised description must not be guessed
    at from its spelling.
    """
    key = " ".join((description or "").lower().split())
    return CLINVAR_CLASSIFICATIONS.get(key)


def _locations(record: dict) -> list[dict]:
    out: list[dict] = []
    for vset in record.get("variation_set") or []:
        out.extend(vset.get("variation_loc") or [])
    return out


def _same_chrom(a: str, b: str) -> bool:
    return str(a).removeprefix("chr").removeprefix("Chr") == str(b).removeprefix("chr").removeprefix("Chr")


def matching_clinvar_records(variant: Variant, records: list[dict]) -> list[dict]:
    """Only the records that are about THIS variant.

    A ClinVar gene search returns every submission for the gene. Applying
    `records[0]` gave a benign variant its neighbour's pathogenic verdict (#61).
    A record matches only when chromosome, position and one of the variant's ALT
    alleles all agree, on the expected assembly.
    """
    alts = {str(a).upper() for a in (variant.alt or [])}
    matched: list[dict] = []

    for record in records:
        for loc in _locations(record):
            if loc.get("assembly_name") and loc["assembly_name"] != CLINVAR_ASSEMBLY:
                continue
            if not _same_chrom(loc.get("chr", ""), variant.chrom):
                continue
            try:
                if int(loc.get("start", -1)) != int(variant.pos):
                    continue
            except (TypeError, ValueError):
                continue
            ref = str(loc.get("ref", "")).upper()
            alt = str(loc.get("alt", "")).upper()
            if ref and ref != str(variant.ref).upper():
                continue
            if alt and alts and alt not in alts:
                continue
            matched.append(record)
            break

    return matched


def significance_for_variant(
    variant: Variant, records: list[dict]
) -> ClinicalSignificance | None:
    """This variant's own ClinVar classification, or None.

    None means "ClinVar has nothing for this variant", which is a real and useful
    answer. It is never a neighbour's verdict.
    """
    for record in matching_clinvar_records(variant, records):
        described = (
            record.get("germline_classification", {}).get("description")
            or record.get("clinical_significance", {}).get("description")
            or ""
        )
        classification = classify_significance(described)
        if classification is not None:
            return classification
    return None


@router.post("/annotate", response_model=VariantAnnotationResponse)
def annotate_variants(req: VariantAnnotationRequest) -> VariantAnnotationResponse:
    annotated_count = 0
    clinvar_errors: list[str] = []
    for v in req.variants:
        changed = False
        if req.annotate_gnomad:
            # Query each ALT allele separately — gnomAD expects one allele per lookup
            for alt_allele in v.alt:
                variant_id = f"{v.chrom}-{v.pos}-{v.ref}-{alt_allele}"
                try:
                    af = _gnomad.get_af(variant_id)
                    if af is not None:
                        v.gnomad_af = af
                        v.allele_frequency = af
                        changed = True
                        break  # use AF from first matched allele
                except Exception:
                    pass
        if req.annotate_clinvar and v.gene:
            try:
                # Still a gene search — ClinVar has no single position endpoint in
                # this client — but the result is now FILTERED to records that are
                # about this variant, and max_results is high enough that the
                # right record is not cut off by an arbitrary 5. A gene with more
                # submissions than this simply yields no match, which is reported
                # as "not found" rather than as a neighbour's verdict (#61).
                records = _ncbi.search_clinvar(v.gene, max_results=CLINVAR_MAX_RECORDS)
                classification = significance_for_variant(v, records)
                if classification is not None:
                    v.clinical_significance = classification
                    changed = True
            except Exception as e:  # noqa: BLE001
                # Was `pass`, so a broken ClinVar lookup was indistinguishable
                # from "no annotation found".
                clinvar_errors.append(f"{v.chrom}:{v.pos} {v.gene}: {type(e).__name__}: {e}")
        if changed:
            annotated_count += 1

    if clinvar_errors:
        # A lookup that failed is not a variant with no annotation. Say so rather
        # than returning a quietly-incomplete result (#61).
        raise HTTPException(
            502,
            "ClinVar lookups failed for "
            f"{len(clinvar_errors)} of {len(req.variants)} variants: "
            + "; ".join(clinvar_errors[:5]),
        )

    return VariantAnnotationResponse(
        variants=req.variants,
        total=len(req.variants),
        annotated=annotated_count,
    )


@router.post("/stats")
def variant_stats(variants: list[Variant]) -> dict:
    """Compute summary statistics for a set of variants."""
    if not variants:
        return {}
    from collections import Counter
    types = Counter(v.variant_type.value for v in variants)
    chroms = Counter(v.chrom for v in variants)
    afs = [v.allele_frequency for v in variants if v.allele_frequency is not None]
    return {
        "total": len(variants),
        "by_type": dict(types),
        "by_chromosome": dict(chroms.most_common(10)),
        "allele_frequency": {
            "mean": round(sum(afs) / len(afs), 6) if afs else None,
            "min": round(min(afs), 6) if afs else None,
            "max": round(max(afs), 6) if afs else None,
        },
        "with_gene": sum(1 for v in variants if v.gene),
        "with_clinical_sig": sum(1 for v in variants if v.clinical_significance),
    }


@router.get("/gnomad/{variant_id}")
def gnomad_variant(variant_id: str) -> dict:
    """Look up a variant in gnomAD (format: CHROM-POS-REF-ALT)."""
    try:
        data = _gnomad.variant(variant_id)
        return data
    except Exception as e:
        raise HTTPException(502, f"gnomAD lookup failed: {e}")
