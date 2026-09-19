"""OR-10, OR-11 — the two spellings of one molecule: RNA, and pasted layout.

Neither test needs a new oracle. Both assert that two strings denoting the SAME
molecule get the same answer, and the answer for the reference spelling has
already been pinned against an independent ORF finder by OR-7.
"""

from __future__ import annotations

import pytest

RNA_CASES = [
    pytest.param("puc19_L09137.2", 2686, 13, id="puc19"),
    pytest.param("NG_007400.1", 24544, 171, id="NG_007400.1"),
]


@pytest.mark.parametrize("name,length,expected_orfs", RNA_CASES)
def test_or_10_rna_input_finds_the_same_orfs_as_its_dna_spelling(
    name, length, expected_orfs, corpus, find_orfs
) -> None:
    """``AUG`` is ``ATG``. An RNA transcript is not a molecule with no ORFs.

    The app models RNA as a first-class molecule type
    (``backend.models.sequence.MoleculeType.RNA``, and ``/sequences/{id}/edit``
    ``transcribe`` produces one), and Biopython's tables are registered over
    both alphabets (``unambiguous_dna_by_id`` and ``unambiguous_rna_by_id``).
    Returning an empty result for a transcript is neither an answer nor a
    refusal.
    """
    dna = corpus[name]
    assert len(dna) == length
    rna = dna.replace("T", "U")

    reference = find_orfs(dna, min_length_aa=30)
    assert reference.total_found == expected_orfs

    result = find_orfs(rna, min_length_aa=30)
    assert result.total_found > 0, (
        f"{name} spelled as RNA reports 0 ORFs; the same molecule spelled as DNA "
        f"reports {expected_orfs}. No error was raised."
    )
    assert {(o.frame, o.start, o.end, o.length_aa, o.protein) for o in result.orfs} == {
        (o.frame, o.start, o.end, o.length_aa, o.protein) for o in reference.orfs
    }


def test_or_11_crlf_wrapped_sequence_denotes_the_same_molecule(col1a1, find_orfs) -> None:
    """NG_007400.1 wrapped at 60 columns with CRLF is still 24,544 bases.

    ``backend/sequence_text.py`` is the app's own definition of which characters
    are the molecule and which are layout — it exists because two endpoints
    disagreeing about that gave one paste two sets of coordinates (#86).
    ``find_orfs`` strips only ``' '`` and ``'\\n'``
    (backend/routers/analysis.py:499), so it is not using it.
    """
    from backend.sequence_text import normalize_template

    wrapped = "\r\n".join(col1a1[i : i + 60] for i in range(0, len(col1a1), 60))
    assert normalize_template(wrapped) == col1a1

    reference = find_orfs(col1a1, min_length_aa=30)
    result = find_orfs(wrapped, min_length_aa=30)

    assert result.sequence_length == 24544, (
        f"sequence_length {result.sequence_length} is the length of no molecule: "
        f"the 409 carriage returns are still in the template"
    )
    assert {(o.frame, o.start, o.end, o.protein) for o in result.orfs} == {
        (o.frame, o.start, o.end, o.protein) for o in reference.orfs
    }
