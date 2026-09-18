"""Structure lookup/prediction + putative interaction points.

Lookup order for /predict: (1) real experimental structure from RCSB PDB,
(2) precomputed AlphaFold DB prediction by UniProt accession, (3) ESMFold
live-fold from the raw sequence as a last resort for sequences with no
UniProt match (novel/synthetic/in-app-edited sequences).
"""

from fastapi import APIRouter, HTTPException

from ..models.structure import (
    StructureSource, StructureRequest, StructureResult,
    InteractionPoint, InteractionPointsRequest, InteractionPointsResult,
)
from ..services import PDBClient, AlphaFoldDBClient, ESMFoldClient
from ..services.interaction_points import compute_interaction_points

router = APIRouter(prefix="/structure", tags=["structure"])

_pdb = PDBClient()
_alphafold_db = AlphaFoldDBClient()
_esmfold = ESMFoldClient()


@router.post("/predict", response_model=StructureResult)
def get_or_predict_structure(req: StructureRequest) -> StructureResult:
    warnings: list[str] = []

    if not req.skip_experimental_lookup and (req.gene_symbol or req.uniprot_id):
        try:
            ids = (
                _pdb.search_by_uniprot(req.uniprot_id) if req.uniprot_id
                else _pdb.search_by_gene(req.gene_symbol)
            )
            if ids:
                chosen = ids[0]
                pdb_text = _pdb.download_pdb(chosen)
                # Say how this entry was picked. RCSB returns result_set in its
                # own relevance order, so ids[0] is "best text match", which is
                # NOT the same as the best structure for the user's sequence:
                # there is no organism filter and no sequence-identity check
                # (#66.3). Claiming otherwise, or saying nothing, leaves the user
                # unable to tell a good match from an arbitrary one.
                basis = (
                    f"Chosen as the top RCSB relevance hit for "
                    f"{'UniProt ' + req.uniprot_id if req.uniprot_id else 'gene ' + str(req.gene_symbol)}"
                )
                if len(ids) > 1:
                    basis += f", out of {len(ids)} entries ({', '.join(ids[1:6])}…)"
                return StructureResult(
                    source=StructureSource.EXPERIMENTAL_PDB,
                    pdb_id=chosen,
                    pdb_text=pdb_text,
                    sequence_length=len(req.sequence),
                    pdb_candidates=ids,
                    confidence_note=(
                        "Experimental structure from RCSB PDB — no per-residue confidence score. "
                        f"{basis}. Not ranked by sequence identity to your sequence, and not "
                        "filtered by organism, so check the entry matches what you mean."
                    ),
                )
        except Exception as e:
            warnings.append(f"PDB lookup failed, trying AlphaFold DB: {e}")

    if req.uniprot_id:
        try:
            prediction = _alphafold_db.get_prediction(req.uniprot_id)
            if prediction is not None:
                return StructureResult(
                    source=StructureSource.PREDICTED_ALPHAFOLD_DB,
                    pdb_id=req.uniprot_id,
                    pdb_text=prediction.pdb_text,
                    sequence_length=len(req.sequence),
                    confidence_note=(
                        "Predicted structure (AlphaFold DB) — B-factor column encodes "
                        "per-residue pLDDT confidence (0-100)."
                        + (
                            f" Overall pLDDT: {prediction.global_metric_value:.1f}."
                            if prediction.global_metric_value is not None else ""
                        )
                    ),
                    warnings=warnings,
                )
            warnings.append("No AlphaFold DB entry for this accession, falling back to ESMFold.")
        except Exception as e:
            warnings.append(f"AlphaFold DB lookup failed, falling back to ESMFold: {e}")

    try:
        pdb_text = _esmfold.predict(req.sequence)
    except ValueError as e:
        raise HTTPException(400, str(e))
    except Exception as e:
        raise HTTPException(502, f"ESMFold prediction failed: {e}")

    return StructureResult(
        source=StructureSource.PREDICTED_ESMFOLD,
        pdb_text=pdb_text,
        sequence_length=len(req.sequence),
        confidence_note="Predicted structure (ESMFold) — B-factor column encodes per-residue pLDDT confidence (0-100).",
        warnings=warnings,
    )


@router.post("/interaction_points", response_model=InteractionPointsResult)
def interaction_points(req: InteractionPointsRequest) -> InteractionPointsResult:
    try:
        points = compute_interaction_points(req.pdb_text)
    except Exception as e:
        raise HTTPException(400, f"Could not parse structure: {e}")
    return InteractionPointsResult(points=[InteractionPoint(**p) for p in points])
