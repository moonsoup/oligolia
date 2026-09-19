"""OR-13, OR-14 — the alternative genetic codes, against NCBI's own proteins.

NC_012920.1 is the revised Cambridge Reference Sequence: the human mitochondrial
genome, 16,569 bp, circular, every one of its 13 CDS features declaring
``/transl_table=2``. Each also carries a ``/translation`` qualifier, which is
NCBI's own protein — a recorded literal that nothing in this repository
produced, and the strongest oracle available to this round.
"""

from __future__ import annotations

import pytest
from orf_oracle import as_initiator_met, ncbi_amino_acids, ncbi_start_codons

# The gene names as NC_012920.1 spells them in its own /gene qualifiers.
GENES = [
    "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3",
    "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB",
]

ND5_SPAN = (12336, 14148)
ND5_LENGTH_AA = 603


def _literal(ann) -> str:
    value = ann.qualifiers.get("translation")
    assert value, f"{ann.qualifiers.get('gene')} carries no /translation"
    return str(value)


def test_or_14_all_thirteen_gene_names_are_present(mito_cds) -> None:
    _, by_gene = mito_cds
    assert sorted(by_gene) == sorted(GENES), sorted(by_gene)


@pytest.mark.parametrize("gene", GENES)
def test_or_14_restated_translation_reproduces_ncbis_protein(gene, mito_cds) -> None:
    """Byte for byte, under the /transl_table=2 the feature declares.

    Covers the minus strand (MT-ND6), the six spans that are not a multiple of
    three because the stop codon is completed post-transcriptionally, and the
    initiator-methionine convention: MT-ND2 opens on ATT, which gc.prt's
    ``sncbieaa`` marks as a table-2 START (``M``), and which NCBI therefore
    writes as M rather than as the I that ATT codes for mid-chain.
    """
    from backend.services.annotations import restated_translation

    sequence, by_gene = mito_cds
    ann = by_gene[gene]
    assert str(ann.qualifiers.get("transl_table")) == "2"

    expected = _literal(ann)
    got = restated_translation(ann, sequence)

    assert got is not None, f"{gene}: restated_translation returned None"
    if got != expected and as_initiator_met(got) == expected:
        start = sequence[ann.start : ann.start + 3] if ann.strand.value != -1 else None
        pytest.fail(
            f"{gene}: initiator not rendered as methionine. NCBI's /translation "
            f"begins {expected[:1]!r}, the app returns {got[:1]!r}; the remaining "
            f"{len(expected) - 1} residues agree. Start codon {start!r} is one of "
            f"table 2's start codons per gc.prt sncbieaa "
            f"({sorted(ncbi_start_codons(2))}), and a start codon is written M."
        )
    assert got == expected, (
        f"{gene}: length {len(got)} vs NCBI's {len(expected)}; "
        f"first difference at residue "
        f"{next((i for i, (a, b) in enumerate(zip(got, expected)) if a != b), 'n/a')}"
    )


def test_or_13_nd5_is_one_orf_under_table_2(mito_cds, find_orfs) -> None:
    """MT-ND5 is a 603-residue reading frame — under its own genetic code.

    Table 2 reads TGA as tryptophan; table 1 reads it as a stop. MT-ND5's 1812
    coding bases contain 11 in-frame TGA codons, so a table-1 finder cannot see
    this gene at all: it shreds it into fragments. ``find_orfs`` takes no
    genetic-code argument (backend/routers/analysis.py:497) and translates with
    a hard-coded table-1 ``GENETIC_CODE`` (:382), so there is no way to ask it.
    """
    sequence, by_gene = mito_cds
    nd5 = by_gene["ND5"]
    assert (nd5.start, nd5.end) == ND5_SPAN
    bases = sequence[nd5.start : nd5.end]
    assert len(bases) == 1812

    codons = [bases[i : i + 3] for i in range(0, len(bases), 3)]
    t1, t2 = ncbi_amino_acids(1), ncbi_amino_acids(2)
    misread = [i for i, c in enumerate(codons) if t1[c] != t2[c]]
    assert codons.count("TGA") == 11
    assert len(_literal(nd5)) == ND5_LENGTH_AA

    import inspect

    params = inspect.signature(find_orfs).parameters
    table_param = {"table", "table_id", "genetic_code", "transl_table", "code"} & set(params)
    assert table_param, (
        f"find_orfs{inspect.signature(find_orfs)} offers no genetic-code "
        f"argument, so NC_012920.1's declared /transl_table=2 cannot be "
        f"requested. Under the hard-coded table 1 its MT-ND5 gene "
        f"({ND5_LENGTH_AA} aa) is misread at {len(misread)} codons "
        f"({codons.count('TGA')} of them TGA, read as stop instead of W), and "
        f"the longest ORF reported over its own bases is "
        f"{max((o.length_aa for o in find_orfs(bases, min_length_aa=30).orfs), default=0)} aa."
    )

    kwargs = {next(iter(table_param)): 2}
    result = find_orfs(bases, min_length_aa=30, **kwargs)
    assert ND5_LENGTH_AA in [o.length_aa for o in result.orfs]
