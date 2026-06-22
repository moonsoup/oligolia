"""Variant annotation and analysis endpoints."""

from fastapi import APIRouter, HTTPException
from ..models.variant import Variant, VariantAnnotationRequest, VariantAnnotationResponse, ClinicalSignificance
from ..services import NCBIClient, GnomADClient

router = APIRouter(prefix="/variants", tags=["variants"])

_ncbi = NCBIClient()
_gnomad = GnomADClient()


@router.post("/annotate", response_model=VariantAnnotationResponse)
def annotate_variants(req: VariantAnnotationRequest) -> VariantAnnotationResponse:
    annotated_count = 0
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
                records = _ncbi.search_clinvar(v.gene, max_results=5)
                if records:
                    sig = records[0].get("clinical_significance", {}).get("description", "")
                    for cs in ClinicalSignificance:
                        if cs.value.lower() in sig.lower():
                            v.clinical_significance = cs
                            changed = True
                            break
            except Exception:
                pass
        if changed:
            annotated_count += 1

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
