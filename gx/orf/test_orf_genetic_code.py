"""OR-1, OR-2, OR-6 — the genetic code itself, and where the start set comes from.

Oracle: ``gx/corpus/ncbi/gc.prt``, NCBI's own tables, read as text. Biopython's
``Bio.Data.CodonTable`` is used only after it has been re-rendered back into
gc.prt's own ``ncbieaa``/``sncbieaa`` strings and shown to agree, so neither
copy of the authority is trusted alone.
"""

from __future__ import annotations

import pytest
from orf_oracle import (
    GC_PRT_CODON_ORDER,
    biopython_matches_gc_prt,
    find_orfs_oracle,
    gc_prt_table,
    ncbi_amino_acids,
    ncbi_start_codons,
)


@pytest.mark.parametrize("table_id", [1, 2, 11])
def test_or_1_the_two_authorities_agree_before_either_is_used(table_id: int) -> None:
    """Biopython's tables 1/2/11 re-render into gc.prt's strings exactly."""
    assert biopython_matches_gc_prt(table_id), (
        f"Biopython's table {table_id} does not re-render into gc.prt's "
        f"ncbieaa/sncbieaa; the round's oracle is unsound, not the subject."
    )


def test_or_1_genetic_code_is_ncbi_table_1() -> None:
    """``GENETIC_CODE`` is table 1: all 64 codons, each amino acid from gc.prt."""
    from backend.routers.analysis import GENETIC_CODE

    expected = ncbi_amino_acids(1)
    assert set(GENETIC_CODE) == set(GC_PRT_CODON_ORDER)
    wrong = {c: (GENETIC_CODE[c], expected[c]) for c in expected if GENETIC_CODE[c] != expected[c]}
    assert wrong == {}, f"codons disagreeing with gc.prt table 1 (app, ncbi): {wrong}"
    assert {c for c, a in GENETIC_CODE.items() if a == "*"} == {"TAA", "TAG", "TGA"}


def test_or_2_gc_prt_says_table_2_changes_exactly_four_codons() -> None:
    """The authority's own statement of what table 2 changes, before the subject."""
    t1, t2 = ncbi_amino_acids(1), ncbi_amino_acids(2)
    differing = {c: (t1[c], t2[c]) for c in t1 if t1[c] != t2[c]}
    assert differing == {
        "TGA": ("*", "W"),
        "ATA": ("I", "M"),
        "AGA": ("R", "*"),
        "AGG": ("R", "*"),
    }


def test_or_2_gc_prt_says_table_11_changes_no_amino_acid_only_starts() -> None:
    """Table 11 is table 1's amino acids with a larger start set."""
    t1, t11 = ncbi_amino_acids(1), ncbi_amino_acids(11)
    assert {c for c in t1 if t1[c] != t11[c]} == set()
    assert ncbi_start_codons(11) - ncbi_start_codons(1) == {"GTG", "ATT", "ATC", "ATA"}
    assert ncbi_start_codons(1) - ncbi_start_codons(11) == set()


def test_or_2_subject_table_2_changes_exactly_those_codons(mito_cds) -> None:
    """MT-CYTB under table 2 vs table 1 differs only where the codon does.

    The subject is ``restated_translation``, the one place in the app that reads
    ``/transl_table`` off a feature. Re-declaring MT-CYTB's table as 1 must
    change the protein at exactly the residue positions whose codon is one of
    the four codons gc.prt says the two tables disagree about — no more, no
    fewer, and nowhere else in a 380-residue chain.
    """
    from backend.services.annotations import restated_translation

    sequence, by_gene = mito_cds
    cytb = by_gene["CYTB"]
    bases = sequence[cytb.start : cytb.end]
    codons = [bases[i : i + 3] for i in range(0, len(bases) // 3 * 3, 3)]

    t1, t2 = ncbi_amino_acids(1), ncbi_amino_acids(2)
    swings = {c for c in t1 if t1[c] != t2[c]}
    expected_positions = {i for i, c in enumerate(codons) if c in swings}
    assert expected_positions, "MT-CYTB carries no table-1/table-2 swing codon"

    as_t2 = restated_translation(cytb, sequence)
    as_t1 = restated_translation(cytb.model_copy(update={
        "qualifiers": {**cytb.qualifiers, "transl_table": "1"},
    }), sequence)

    # `restated_translation` passes to_stop=True, so the table-1 rendering ends
    # at MT-CYTB's first codon that table 1 calls a stop — the first in-frame
    # TGA, which table 2 reads as tryptophan and carries straight past.
    t1_stop = min(i for i, c in enumerate(codons) if t1[c] == "*")
    assert len(as_t1) == t1_stop, (
        f"table 1 should stop at MT-CYTB's first table-1 stop codon "
        f"({codons[t1_stop]} at residue {t1_stop}); got {len(as_t1)} residues"
    )

    differing = {i for i, (a, b) in enumerate(zip(as_t1, as_t2)) if a != b}
    should_differ = {i for i in expected_positions if i < t1_stop}
    assert differing == should_differ, (
        f"table 2 changed residues whose codon is not one of {sorted(swings)}: "
        f"{sorted(differing - should_differ)[:10]}; and left unchanged residues "
        f"whose codon is: {sorted(should_differ - differing)[:10]}"
    )


def test_or_2_subject_table_11_equals_table_1(mito_cds) -> None:
    """A /transl_table=11 feature translates identically to /transl_table=1."""
    from backend.services.annotations import restated_translation

    sequence, by_gene = mito_cds
    cytb = by_gene["CYTB"]
    as_t1 = restated_translation(
        cytb.model_copy(update={"qualifiers": {**cytb.qualifiers, "transl_table": "1"}}),
        sequence,
    )
    as_t11 = restated_translation(
        cytb.model_copy(update={"qualifiers": {**cytb.qualifiers, "transl_table": "11"}}),
        sequence,
    )
    assert as_t11 == as_t1


def test_or_6_include_all_starts_uses_the_tables_start_codons(col1a1, find_orfs) -> None:
    """``include_all_starts=True`` must mean 'every start codon of table 1'.

    gc.prt's ``sncbieaa`` for table 1 marks exactly ATG, CTG and TTG. The
    comparison is against an independent finder given that same set, over the
    full 24,544-base NG_007400.1, at min_length_aa=30.
    """
    _, sncbieaa = gc_prt_table(1)
    assert ncbi_start_codons(1) == {"ATG", "CTG", "TTG"}, sncbieaa

    reported = {
        (o.frame, o.start, o.end, o.length_aa, o.start_codon)
        for o in find_orfs(col1a1, min_length_aa=30, include_all_starts=True).orfs
    }
    expected = {
        (o["frame"], o["start"], o["end"], o["length_aa"], o["start_codon"])
        for o in find_orfs_oracle(col1a1, min_length_aa=30, table_id=1)
    }

    ctg_found = {o for o in reported if o[4] == "CTG"}
    ctg_real = {o for o in expected if o[4] == "CTG"}
    non_table_1 = {o for o in reported if o[4] not in ncbi_start_codons(1)}

    assert ctg_found == ctg_real, (
        f"CTG is a table-1 start codon: {len(ctg_real)} CTG-initiated ORFs of "
        f">=30 aa exist in NG_007400.1, {len(ctg_found)} were reported"
    )
    assert non_table_1 == set(), (
        "ORFs opened on codons gc.prt does not mark as table-1 starts: "
        f"{sorted({o[4] for o in non_table_1})}"
    )
    assert reported == expected
