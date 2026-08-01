"""Heuristic putative-interaction-point annotation.

Classifies each residue in a structure by side-chain charge/polarity and combines
that with per-residue relative solvent-accessible surface area (Shrake-Rupley) to
flag surface-exposed charged/polar residues as "putative interaction points" —
residues plausibly available to participate in intermolecular electrostatic/polar
contacts.

This is a classical surface/charge heuristic, NOT a validated binding-site,
docking, or protein-protein-interaction predictor. It has no concept of a second
molecule, geometry of approach, or binding energetics — see the `METHOD_NOTE`
string, which the router surfaces verbatim to the GUI.
"""

from io import StringIO

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

METHOD_NOTE = (
    "Charge/polarity classification (static side-chain lookup table) + Shrake-Rupley "
    "relative solvent-accessible surface area, thresholded. Heuristic surface/charge "
    "annotation — NOT a validated binding-site or docking prediction."
)

CHARGE_CLASS = {
    "ASP": "acidic", "GLU": "acidic",
    "LYS": "basic", "ARG": "basic", "HIS": "basic",
    "SER": "polar", "THR": "polar", "ASN": "polar",
    "GLN": "polar", "TYR": "polar", "CYS": "polar",
}
# Everything else (ALA, VAL, LEU, ILE, PHE, TRP, MET, PRO, GLY) -> "hydrophobic"

# Theoretical maximum solvent-accessible surface area per residue, in Å²
# (Tien et al. 2013, "Maximum Allowed Solvent Accessibilities of Residues in
# Proteins", PLOS ONE) — a small embedded reference table, not a new dependency.
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

SASA_EXPOSURE_THRESHOLD = 0.20  # relative SASA fraction


def compute_interaction_points(pdb_text: str) -> list[dict]:
    """Parse pdb_text, compute per-residue relative SASA, classify, threshold.

    Returns one dict per standard amino-acid residue:
    residue_index, residue_name, chain, classification, relative_sasa,
    is_putative_interaction_point.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", StringIO(pdb_text))

    sr = ShrakeRupley()
    sr.compute(structure, level="R")

    points: list[dict] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                max_asa = MAX_ASA.get(resname)
                if max_asa is None:
                    continue  # skip waters/ligands/non-standard residues
                classification = CHARGE_CLASS.get(resname, "hydrophobic")
                relative_sasa = float(min(residue.sasa / max_asa, 1.0))
                is_interaction_point = bool(
                    classification in ("acidic", "basic", "polar")
                    and relative_sasa >= SASA_EXPOSURE_THRESHOLD
                )
                points.append({
                    "residue_index": int(residue.id[1]),
                    "residue_name": resname,
                    "chain": chain.id,
                    "classification": classification,
                    "relative_sasa": round(relative_sasa, 4),
                    "is_putative_interaction_point": is_interaction_point,
                })
        break  # first model only — structures here are single-model (X-ray/predicted)

    return points
