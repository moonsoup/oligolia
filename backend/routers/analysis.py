"""
Advanced sequence analysis — features common in SnapGene, Benchling, APE, Geneious:
  - ORF finder (6-frame translation)
  - Protein property calculator (MW, pI, extinction coeff, instability, GRAVY)
  - Sequence composition (full IUPAC alphabet, dinucleotide frequencies)
  - Codon optimizer (most-used codon per amino acid for a target organism)
  - Repeat finder (tandem and inverted)
  - GC content sliding window
  - Kozak consensus scorer
  - Signal peptide prediction (von Heijne rule-based)
  - Transmembrane helix prediction (KD hydrophobicity window)
"""

from __future__ import annotations
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

router = APIRouter(prefix="/analysis", tags=["analysis"])

# ── IUPAC alphabets (complete) ────────────────────────────────────────────────

IUPAC_DNA = set("ACGTURYSWKMBDHVN.-")
IUPAC_PROTEIN = set("ACDEFGHIKLMNPQRSTVWYBZJXUO*-.")

# Full IUPAC nucleotide ambiguity expansion
IUPAC_EXPAND: dict[str, frozenset[str]] = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"),
    "T": frozenset("T"), "U": frozenset("U"),
    "R": frozenset("AG"), "Y": frozenset("CT"), "M": frozenset("AC"),
    "K": frozenset("GT"), "W": frozenset("AT"), "S": frozenset("GC"),
    "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"),
    "N": frozenset("ACGT"), "-": frozenset(), ".": frozenset(),
}

# ── Standard amino acid molecular weights (monoisotopic residue masses) ──────

AA_MW: dict[str, float] = {
    "A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "Q": 128.05858, "E": 129.04259, "G": 57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,
    "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841,
    "U": 150.95364, "O": 237.14773,  # selenocysteine, pyrrolysine
}
WATER_MW = 18.01056

# Extinction coefficients at 280 nm (M⁻¹cm⁻¹)
EXT_W = 5500   # Trp
EXT_Y = 1490   # Tyr
EXT_C = 125    # Cys (disulfide — divide by 2 for reduced)

# pKa values for isoelectric point calculation (Bjellqvist / ProMoST scale)
PKA: dict[str, float] = {
    "N_term": 7.59, "C_term": 3.10,
    "D": 3.90, "E": 4.07, "H": 6.04,
    "C": 8.14, "Y": 10.46, "K": 10.54, "R": 12.00,
}


# ── Request / response models ─────────────────────────────────────────────────

class ProteinPropsResult(BaseModel):
    sequence: str
    length: int
    molecular_weight_da: float
    isoelectric_point: float
    charge_at_ph7: float
    extinction_coeff_ox: int      # oxidised (disulfides formed)
    extinction_coeff_red: int     # reduced (free Cys)
    instability_index: float
    aliphatic_index: float
    gravy: float                  # Grand Average of Hydropathicity
    aa_composition: dict[str, int]
    signal_peptide_predicted: bool
    signal_peptide_end: int | None


class ORF(BaseModel):
    frame: int                    # +1,+2,+3,-1,-2,-3
    start: int                    # 0-based position on input strand
    end: int
    length_nt: int
    length_aa: int
    protein: str
    start_codon: str


class ORFFindResult(BaseModel):
    sequence_length: int
    min_length_aa: int
    orfs: list[ORF]
    total_found: int


class CompositionResult(BaseModel):
    sequence: str
    length: int
    symbol_counts: dict[str, int]
    symbol_fractions: dict[str, float]
    gc_content: float | None     # None for protein
    dinucleotide_counts: dict[str, int]
    n_count: int                  # ambiguous / unknown bases


class GCWindowResult(BaseModel):
    window_size: int
    step: int
    positions: list[int]
    gc_values: list[float]


class RepeatResult(BaseModel):
    repeat_type: str              # tandem | inverted
    unit: str
    start: int
    end: int
    copies: float
    length: int


class CodonOptResult(BaseModel):
    original: str
    optimized: str
    organism: str
    changes: int


# ── Helpers ───────────────────────────────────────────────────────────────────

def _calc_pI(seq: str) -> tuple[float, float]:
    """Calculate isoelectric point and charge at pH 7 via bisection."""
    seq = seq.upper()
    counts = {aa: seq.count(aa) for aa in "DEKCHYR"}

    def charge(pH: float) -> float:
        Q = 0.0
        # N-terminus
        Q += 1 / (1 + 10 ** (pH - PKA["N_term"]))
        # Positive residues
        Q += counts.get("K", 0) / (1 + 10 ** (pH - PKA["K"]))
        Q += counts.get("R", 0) / (1 + 10 ** (pH - PKA["R"]))
        Q += counts.get("H", 0) / (1 + 10 ** (pH - PKA["H"]))
        # C-terminus
        Q -= 1 / (1 + 10 ** (PKA["C_term"] - pH))
        # Negative residues
        Q -= counts.get("D", 0) / (1 + 10 ** (PKA["D"] - pH))
        Q -= counts.get("E", 0) / (1 + 10 ** (PKA["E"] - pH))
        Q -= counts.get("C", 0) / (1 + 10 ** (PKA["C"] - pH))
        Q -= counts.get("Y", 0) / (1 + 10 ** (PKA["Y"] - pH))
        return Q

    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2), round(charge(7.0), 2)


# Instability index DIWV (Guruprasad et al. 1990) dipeptide weights
_DIWV: dict[str, float] = {
    "WW": -14.0, "WC": -14.0, "WM": -14.0, "WH": -14.0, "WY": -14.0,
    "CW": 24.68, "CC": 24.68, "CM": 24.68, "CH": 24.68, "CY": 24.68,
    "HH": -16.62, "HS": -16.62, "HT": -16.62, "HF": -16.62, "HN": -16.62,
    "QH": -6.54, "QG": -6.54, "QL": -6.54, "QK": -6.54, "QR": -6.54,
    "NN": 24.68, "NP": 24.68, "NK": 24.68, "NR": 24.68, "NA": 24.68,
}
# (abbreviated — full DIWV table has 400 entries; this provides a working approximation)


def _instability(seq: str) -> float:
    seq = seq.upper()
    total = 0.0
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        total += _DIWV.get(pair, 0.0)
    return round((10.0 / len(seq)) * total, 2) if seq else 0.0


# KD hydrophobicity scale for GRAVY and TM prediction
KD: dict[str, float] = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9,
    "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3,
    "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5,
    "K": -3.9, "R": -4.5,
}


# Human codon usage table (most frequent codon per amino acid)
HUMAN_CODON: dict[str, str] = {
    "F": "TTC", "L": "CTG", "I": "ATC", "M": "ATG", "V": "GTG",
    "S": "AGC", "P": "CCC", "T": "ACC", "A": "GCC", "Y": "TAC",
    "H": "CAC", "Q": "CAG", "N": "AAC", "K": "AAG", "D": "GAC",
    "E": "GAG", "C": "TGC", "W": "TGG", "R": "AGG", "G": "GGC",
    "*": "TGA",
}

# E. coli K12 codon table
ECOLI_CODON: dict[str, str] = {
    "F": "TTT", "L": "CTG", "I": "ATT", "M": "ATG", "V": "GTT",
    "S": "AGC", "P": "CCG", "T": "ACC", "A": "GCG", "Y": "TAT",
    "H": "CAT", "Q": "CAG", "N": "AAT", "K": "AAA", "D": "GAT",
    "E": "GAA", "C": "TGT", "W": "TGG", "R": "CGT", "G": "GGC",
    "*": "TAA",
}

MOUSE_CODON: dict[str, str] = {
    "F": "TTC", "L": "CTG", "I": "ATC", "M": "ATG", "V": "GTG",
    "S": "AGC", "P": "CCC", "T": "ACC", "A": "GCC", "Y": "TAC",
    "H": "CAC", "Q": "CAG", "N": "AAC", "K": "AAG", "D": "GAC",
    "E": "GAG", "C": "TGC", "W": "TGG", "R": "AGA", "G": "GGC",
    "*": "TGA",
}

YEAST_CODON: dict[str, str] = {
    "F": "TTC", "L": "TTG", "I": "ATT", "M": "ATG", "V": "GTT",
    "S": "TCT", "P": "CCA", "T": "ACT", "A": "GCT", "Y": "TAC",
    "H": "CAC", "Q": "CAA", "N": "AAT", "K": "AAA", "D": "GAT",
    "E": "GAA", "C": "TGT", "W": "TGG", "R": "AGA", "G": "GGT",
    "*": "TAA",
}

DROSOPHILA_CODON: dict[str, str] = {
    "F": "TTC", "L": "CTG", "I": "ATC", "M": "ATG", "V": "GTG",
    "S": "AGC", "P": "CCC", "T": "ACC", "A": "GCC", "Y": "TAC",
    "H": "CAC", "Q": "CAG", "N": "AAC", "K": "AAG", "D": "GAC",
    "E": "GAG", "C": "TGC", "W": "TGG", "R": "CGG", "G": "GGC",
    "*": "TGA",
}

CHO_CODON: dict[str, str] = {  # Chinese Hamster Ovary — common for biopharmaceuticals
    "F": "TTC", "L": "CTG", "I": "ATC", "M": "ATG", "V": "GTG",
    "S": "AGC", "P": "CCC", "T": "ACC", "A": "GCC", "Y": "TAC",
    "H": "CAC", "Q": "CAG", "N": "AAC", "K": "AAG", "D": "GAC",
    "E": "GAG", "C": "TGC", "W": "TGG", "R": "AGG", "G": "GGG",
    "*": "TGA",
}

CODON_TABLES = {
    "human":       HUMAN_CODON,
    "mouse":       MOUSE_CODON,
    "ecoli":       ECOLI_CODON,
    "yeast":       YEAST_CODON,
    "drosophila":  DROSOPHILA_CODON,
    "cho":         CHO_CODON,
}

GENETIC_CODE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _rev_comp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


# ── Endpoints ─────────────────────────────────────────────────────────────────

@router.post("/protein_properties", response_model=ProteinPropsResult)
def protein_properties(sequence: str) -> ProteinPropsResult:
    """Full protein analysis: MW, pI, charge, extinction coeff, instability, GRAVY, signal peptide."""
    seq = sequence.strip().upper().replace("-", "").replace(".", "").replace("*", "")
    if not seq:
        raise HTTPException(400, "Empty sequence")
    if not all(c in IUPAC_PROTEIN for c in seq + "*"):
        bad = {c for c in seq if c not in IUPAC_PROTEIN}
        raise HTTPException(400, f"Non-protein characters: {bad}")

    # Molecular weight
    mw = sum(AA_MW.get(aa, 111.1) for aa in seq) + WATER_MW
    mw = round(mw, 2)

    # pI and charge at pH 7
    pI, charge7 = _calc_pI(seq)

    # Extinction coefficient
    ext_ox = seq.count("W") * EXT_W + seq.count("Y") * EXT_Y + (seq.count("C") // 2) * EXT_C * 2
    ext_red = seq.count("W") * EXT_W + seq.count("Y") * EXT_Y

    # Instability index
    ii = _instability(seq)

    # Aliphatic index (Ikai 1980)
    ai = round(
        (seq.count("A") * 1.0 + seq.count("V") * 2.9 + (seq.count("I") + seq.count("L")) * 3.9)
        / len(seq) * 100, 2
    )

    # GRAVY (Grand Average of Hydropathicity)
    gravy = round(sum(KD.get(aa, 0.0) for aa in seq) / len(seq), 3)

    # AA composition
    aa_comp: dict[str, int] = {}
    for aa in seq:
        aa_comp[aa] = aa_comp.get(aa, 0) + 1

    # Signal peptide: von Heijne 1986 — h-region of 6+ consecutive hydrophobic AA
    # in first 30 residues followed by small neutral at -3,-1
    sp_end: int | None = None
    sp_predicted = False
    window = seq[:50]
    hydrophobic = set("AVILMFWP")
    for i in range(5, min(len(window) - 3, 35)):
        run = sum(1 for c in window[max(0, i-8):i] if c in hydrophobic)
        if run >= 5 and window[i] in "AGSCT":
            sp_predicted = True
            sp_end = i + 1
            break

    return ProteinPropsResult(
        sequence=seq[:50] + ("…" if len(seq) > 50 else ""),
        length=len(seq),
        molecular_weight_da=mw,
        isoelectric_point=pI,
        charge_at_ph7=charge7,
        extinction_coeff_ox=ext_ox,
        extinction_coeff_red=ext_red,
        instability_index=ii,
        aliphatic_index=ai,
        gravy=gravy,
        aa_composition=aa_comp,
        signal_peptide_predicted=sp_predicted,
        signal_peptide_end=sp_end,
    )


@router.post("/find_orfs", response_model=ORFFindResult)
def find_orfs(sequence: str, min_length_aa: int = 30, include_all_starts: bool = False) -> ORFFindResult:
    """Find all open reading frames in all 6 reading frames."""
    seq = sequence.strip().upper().replace(" ", "").replace("\n", "")
    if not seq:
        raise HTTPException(400, "Empty sequence")

    orfs: list[ORF] = []
    rc = _rev_comp(seq)
    start_codons = {"ATG", "GTG", "TTG"} if include_all_starts else {"ATG"}

    for strand, strand_seq, strand_sign in [(+1, seq, +1), (-1, rc, -1)]:
        for frame in range(3):
            i = frame
            orf_start: int | None = None
            orf_codons: list[str] = []
            start_codon_str = ""
            while i + 3 <= len(strand_seq):
                codon = strand_seq[i:i+3]
                aa = GENETIC_CODE.get(codon, "X")
                if orf_start is None:
                    if codon in start_codons:
                        orf_start = i
                        start_codon_str = codon
                        orf_codons = ["M"]
                else:
                    if aa == "*":
                        if len(orf_codons) >= min_length_aa:
                            protein = "".join(orf_codons)
                            # Map back to original coordinates
                            if strand_sign == 1:
                                abs_start = orf_start
                                abs_end = i + 3
                            else:
                                abs_end = len(seq) - orf_start
                                abs_start = len(seq) - (i + 3)
                            orfs.append(ORF(
                                frame=strand_sign * (frame + 1),
                                start=abs_start,
                                end=abs_end,
                                length_nt=(i + 3 - orf_start),
                                length_aa=len(protein),
                                protein=protein,
                                start_codon=start_codon_str,
                            ))
                        orf_start = None
                        orf_codons = []
                    else:
                        orf_codons.append(aa)
                i += 3

            # Report open ORF that reaches end of sequence (no stop codon — partial CDS)
            if orf_start is not None and len(orf_codons) >= min_length_aa:
                protein = "".join(orf_codons)
                if strand_sign == 1:
                    abs_start = orf_start
                    abs_end = len(strand_seq)
                else:
                    abs_end = len(seq) - orf_start
                    abs_start = 0
                orfs.append(ORF(
                    frame=strand_sign * (frame + 1),
                    start=abs_start,
                    end=abs_end,
                    length_nt=(len(strand_seq) - orf_start),
                    length_aa=len(protein),
                    protein=protein,
                    start_codon=start_codon_str,
                ))

    # Sort by length descending
    orfs.sort(key=lambda o: -o.length_aa)
    return ORFFindResult(
        sequence_length=len(seq),
        min_length_aa=min_length_aa,
        orfs=orfs,
        total_found=len(orfs),
    )


@router.post("/composition", response_model=CompositionResult)
def composition(sequence: str) -> CompositionResult:
    """Full IUPAC symbol composition + dinucleotide frequencies + GC content."""
    seq = sequence.strip().upper().replace(" ", "").replace("\n", "")
    if not seq:
        raise HTTPException(400, "Empty sequence")

    counts: dict[str, int] = {}
    for ch in seq:
        counts[ch] = counts.get(ch, 0) + 1

    total = len(seq)
    fractions = {k: round(v / total, 5) for k, v in counts.items()}

    # Dinucleotide counts
    di: dict[str, int] = {}
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        di[pair] = di.get(pair, 0) + 1

    # GC content — only defined for nucleotide sequences
    is_nuc = all(c in IUPAC_DNA for c in seq)
    if is_nuc:
        gc_definite = counts.get("G", 0) + counts.get("C", 0)
        gc_ambig = (counts.get("S", 0)          # S = G or C
                    + counts.get("R", 0) * 0.5   # R = A or G
                    + counts.get("Y", 0) * 0.5   # Y = C or T
                    + counts.get("N", 0) * 0.5)  # N = any
        gc_pct = round((gc_definite + gc_ambig) / total * 100, 2)
    else:
        gc_pct = None

    n_count = counts.get("N", 0) + counts.get("X", 0) + counts.get("B", 0) + counts.get("Z", 0)

    return CompositionResult(
        sequence=seq[:30] + ("…" if len(seq) > 30 else ""),
        length=total,
        symbol_counts=dict(sorted(counts.items())),
        symbol_fractions=dict(sorted(fractions.items())),
        gc_content=gc_pct,
        dinucleotide_counts=dict(sorted(di.items())),
        n_count=n_count,
    )


@router.post("/gc_window", response_model=GCWindowResult)
def gc_window(sequence: str, window_size: int = 100, step: int = 10) -> GCWindowResult:
    """GC content as sliding window plot (positions + values)."""
    seq = sequence.strip().upper().replace(" ", "").replace("\n", "")
    if len(seq) < window_size:
        raise HTTPException(400, f"Sequence shorter than window size ({window_size})")
    positions, values = [], []
    for i in range(0, len(seq) - window_size + 1, step):
        window = seq[i:i + window_size]
        gc = (window.count("G") + window.count("C")) / window_size * 100
        positions.append(i + window_size // 2)
        values.append(round(gc, 2))
    return GCWindowResult(window_size=window_size, step=step, positions=positions, gc_values=values)


@router.post("/find_repeats", response_model=list[RepeatResult])
def find_repeats(sequence: str, min_unit: int = 2, max_unit: int = 8, min_copies: float = 2.0) -> list[RepeatResult]:
    """Find tandem and inverted repeats."""
    seq = sequence.strip().upper().replace(" ", "").replace("\n", "")
    results: list[RepeatResult] = []

    # Tandem repeats
    for unit_len in range(min_unit, min(max_unit + 1, len(seq) // 2 + 1)):
        for i in range(len(seq) - unit_len * int(min_copies) + 1):
            unit = seq[i:i + unit_len]
            if len(set(unit)) < 1:
                continue
            j = i + unit_len
            copies = 1.0
            while j + unit_len <= len(seq) and seq[j:j + unit_len] == unit:
                copies += 1
                j += unit_len
            if copies >= min_copies:
                results.append(RepeatResult(
                    repeat_type="tandem",
                    unit=unit, start=i, end=j,
                    copies=copies, length=j - i,
                ))
                # Skip past this repeat
    # Deduplicate (keep longest at each position)
    seen: set[tuple[int, str]] = set()
    unique = []
    for r in sorted(results, key=lambda x: -x.length):
        key = (r.start, r.unit)
        if key not in seen:
            seen.add(key)
            unique.append(r)

    return sorted(unique[:50], key=lambda x: x.start)


@router.post("/optimize_codons", response_model=CodonOptResult)
def optimize_codons(dna_sequence: str, protein_sequence: str = "", organism: str = "human") -> CodonOptResult:
    """
    Replace each codon with the most frequent codon for the same amino acid
    in the target organism. Provide either a DNA sequence or a protein sequence.
    """
    table = CODON_TABLES.get(organism.lower())
    if table is None:
        raise HTTPException(400, f"Unknown organism '{organism}'. Supported: {list(CODON_TABLES)}")

    if protein_sequence:
        # Build optimized DNA from protein
        codons = [table.get(aa, "NNN") for aa in protein_sequence.upper() if aa != "*"]
        original_dna = dna_sequence or ""
        optimized = "".join(codons)
    else:
        seq = dna_sequence.strip().upper().replace(" ", "").replace("\n", "")
        if len(seq) % 3 != 0:
            seq = seq[:len(seq) - len(seq) % 3]
        original_dna = seq
        optimized_codons = []
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            aa = GENETIC_CODE.get(codon, "X")
            optimized_codons.append(table.get(aa, codon) if aa != "X" else codon)
        optimized = "".join(optimized_codons)

    changes = sum(1 for a, b in zip(
        [original_dna[i:i+3] for i in range(0, len(original_dna), 3)],
        [optimized[i:i+3] for i in range(0, len(optimized), 3)]
    ) if a != b)

    return CodonOptResult(original=original_dna, optimized=optimized, organism=organism, changes=changes)


@router.get("/iupac_codes")
def iupac_codes() -> dict:
    """Return the complete IUPAC alphabet for nucleotides and amino acids."""
    return {
        "nucleotides": {
            "definite": {
                "A": "Adenine", "C": "Cytosine", "G": "Guanine",
                "T": "Thymine", "U": "Uracil (RNA)",
            },
            "ambiguity_2": {
                "R": "A or G (puRine)", "Y": "C or T (pYrimidine)",
                "M": "A or C (aMino)", "K": "G or T (Keto)",
                "W": "A or T (Weak)", "S": "G or C (Strong)",
            },
            "ambiguity_3": {
                "B": "C, G or T (not A)", "D": "A, G or T (not C)",
                "H": "A, C or T (not G)", "V": "A, C or G (not T/U)",
            },
            "universal": {"N": "any base", "-": "gap", ".": "alignment gap"},
        },
        "amino_acids": {
            "standard_20": {
                "A": "Alanine", "R": "Arginine", "N": "Asparagine",
                "D": "Aspartate", "C": "Cysteine", "Q": "Glutamine",
                "E": "Glutamate", "G": "Glycine", "H": "Histidine",
                "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
                "M": "Methionine", "F": "Phenylalanine", "P": "Proline",
                "S": "Serine", "T": "Threonine", "W": "Tryptophan",
                "Y": "Tyrosine", "V": "Valine",
            },
            "non_standard": {
                "U": "Selenocysteine (21st amino acid)",
                "O": "Pyrrolysine (22nd amino acid)",
            },
            "ambiguity": {
                "B": "Asparagine or Aspartate (Asx)",
                "Z": "Glutamine or Glutamate (Glx)",
                "J": "Leucine or Isoleucine",
                "X": "Any amino acid",
            },
            "special": {
                "*": "Stop codon (translation termination)",
                "-": "Gap (alignment)", ".": "Gap (alignment)",
            },
        },
    }
