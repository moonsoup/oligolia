"""OR-15 — the one surface that hands a user a translated sequence.

``POST /sequences/{id}/edit`` with ``operation='translate'``
(backend/routers/sequences.py:115) calls ``Seq(seq_str).translate(to_stop=True)``
on the whole stored sequence. The GUI duplicates the same call at
gui/panels/sequence_panel.py:1077.

This is where the backend suite's ``BiopythonWarning: Partial codon`` comes
from. The question this round asks is whether the app handles that warning or
passes it through, and the answer has to be read off the app's own result: a
desktop user never sees a Python warning.

The DNA used is real: the first 17 bases of pUC19 as pinned, not a constructed
string. Any 17 bases would do — the point is the length.
"""

from __future__ import annotations

from Bio import BiopythonWarning
from Bio.Seq import Seq
from orf_oracle import biopython_translation, partial_tail


def test_or_15_a_ragged_input_is_refused_or_reported(puc19, translate_op) -> None:
    """17 bases cannot become 5 residues without 2 bases being dropped."""
    bases = puc19[:17]
    assert len(bases) == 17

    result, caught = translate_op(bases)

    biopython_warned = [
        str(w.message) for w in caught if issubclass(w.category, BiopythonWarning)
    ]
    expected_protein = biopython_translation(bases, 1)
    assert result.result_seq == expected_protein
    assert len(expected_protein) == 5 and partial_tail(bases) == bases[15:]

    reported = (result.message or "").lower()
    mentions = any(
        token in reported
        for token in ("partial", "trim", "dropped", "incomplete", "not a multiple", "remainder")
    )
    assert not biopython_warned or mentions, (
        f"translate raised BiopythonWarning({biopython_warned[0]!r}) and reported "
        f"{result.message!r} — a message that states 17 nt became 5 aa without "
        f"saying that {len(partial_tail(bases))} base(s) ({partial_tail(bases)!r}) "
        f"were silently discarded. A GUI user sees only the message."
    )


def test_or_15_a_clean_input_raises_no_warning(puc19, translate_op) -> None:
    """The control: 15 bases is five whole codons and must be quiet."""
    bases = puc19[:15]
    result, caught = translate_op(bases)
    assert [str(w.message) for w in caught if issubclass(w.category, BiopythonWarning)] == []
    assert result.result_seq == biopython_translation(bases, 1)


def test_or_15_translation_matches_biopython_on_whole_codons(puc19, translate_op) -> None:
    """Over a real 2685-base whole-codon prefix of pUC19, table 1, to_stop."""
    bases = puc19[:2685]
    assert len(bases) % 3 == 0
    result, caught = translate_op(bases)

    full = str(Seq(bases).translate())
    expected = full.split("*")[0]
    assert result.result_seq == expected, (
        f"translate returned {len(result.result_seq)} residues; "
        f"Bio.Seq.translate(to_stop=True) on the same bases gives {len(expected)}"
    )
    assert [str(w.message) for w in caught if issubclass(w.category, BiopythonWarning)] == []
