"""OR-9 — partial trailing codons, on real bases NCBI itself calls incomplete.

The two cases are the rCRS CDS features whose spans are not a multiple of three
because the stop codon is completed post-transcriptionally by polyadenylation:
MT-ND4 ``[10759:12137)``, 1378 nt, ``/transl_except=(pos:12137,aa:TERM)``, and
MT-CYTB ``[14746:15887)``, 1141 nt, ``/transl_except=(pos:15887,aa:TERM)``.
Handing those bases to the ORF finder is the real form of the question, not a
constructed one: a reading frame that reaches the end of the molecule without a
stop codon is exactly what a partial CDS is.
"""

from __future__ import annotations

import warnings

import pytest
from Bio import BiopythonWarning
from Bio.Seq import Seq
from orf_oracle import ncbi_amino_acids, orf_bases, partial_tail

STOPS = {c for c, a in ncbi_amino_acids(1).items() if a == "*"}

RAGGED_CDS = [
    pytest.param("ND4", 1378, id="ND4"),
    pytest.param("CYTB", 1141, id="CYTB"),
]


@pytest.fixture
def ragged_bases(mito_cds):
    def _bases(gene: str) -> str:
        sequence, by_gene = mito_cds
        ann = by_gene[gene]
        return sequence[ann.start : ann.end]

    return _bases


@pytest.mark.parametrize("gene,length_nt", RAGGED_CDS)
def test_or_9_re_extracting_an_orf_does_not_raise_partial_codon(
    gene, length_nt, ragged_bases, find_orfs
) -> None:
    """No reported interval, sliced back out, is a partial codon to Biopython.

    ``Bio.Seq.translate`` warns 'Partial codon, len(sequence) not a multiple of
    three. Explicitly trim the sequence or add trailing N before translation.
    This may become an error in future.' Any consumer of an ORF interval — the
    app's own Assembly and workflow paths included — does exactly this slice.
    """
    bases = ragged_bases(gene)
    assert len(bases) == length_nt and length_nt % 3 != 0

    offenders = []
    for orf in find_orfs(bases, min_length_aa=30).orfs:
        extracted = orf_bases(bases, orf.start, orf.end, orf.frame)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            Seq(extracted).translate()
        if any(issubclass(w.category, BiopythonWarning) for w in caught):
            offenders.append(
                (orf.frame, orf.start, orf.end, orf.length_nt, orf.length_aa,
                 partial_tail(extracted))
            )

    assert offenders == [], (
        f"{gene}: reported ORF interval(s) whose bases are not a whole number of "
        f"codons (frame, start, end, length_nt, length_aa, untranslated tail): "
        f"{offenders}"
    )


@pytest.mark.parametrize("gene,length_nt", RAGGED_CDS)
def test_or_9_an_unterminated_orf_is_distinguishable_from_a_complete_one(
    gene, length_nt, ragged_bases, find_orfs
) -> None:
    """A run-off ORF must say so; a consumer must not have to re-derive it.

    The subject reports a run-off frame under the same shape as a terminated
    one: ``protein`` is unterminated either way, and ``ORF`` carries no field
    distinguishing them (backend/routers/analysis.py:110). The only difference
    is arithmetic on ``length_nt``, which requires knowing the stop-inclusive
    convention — and the reported bases do not even contain the residues the
    arithmetic would describe.
    """
    bases = ragged_bases(gene)
    result = find_orfs(bases, min_length_aa=30)

    run_off = []
    for orf in result.orfs:
        extracted = orf_bases(bases, orf.start, orf.end, orf.frame)
        codons = [extracted[i : i + 3] for i in range(0, len(extracted) // 3 * 3, 3)]
        if codons and codons[-1] not in STOPS:
            run_off.append(orf)
    assert run_off, f"{gene}: expected at least one frame to run off the end"

    fields = set(type(run_off[0]).model_fields)
    assert fields & {"partial", "is_partial", "truncated", "has_stop", "complete"}, (
        f"{gene}: {len(run_off)} ORF(s) reach the end of the molecule without a "
        f"stop codon, and the ORF model carries no field saying so: "
        f"{sorted(fields)}"
    )
