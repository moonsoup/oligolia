"""Fixtures for the gx translation/ORF round. Oracles live in ``orf_oracle.py``.

Nothing is imported from here by name -- ``from conftest import ...`` is what
forced earlier rounds into separate processes, so this file is fixtures only.

The subject is the app's translation and open-reading-frame code, called
in-process:

* ``backend/routers/analysis.find_orfs`` -- the six-frame ORF finder, its
  hard-coded ``GENETIC_CODE`` dict, ``min_length_aa`` and ``include_all_starts``.
* ``backend/routers/sequences.edit_sequence`` with ``operation="translate"`` --
  the only surface that hands a user a translated frame.
* ``backend/services/annotations.restated_translation`` -- the one place that
  reads ``/transl_table`` and ``/codon_start`` off a feature.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from Bio import SeqIO

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

INPUTS = REPO_ROOT / "gx" / "corpus" / "inputs"


def _genbank_sequence(path: Path) -> str:
    return str(SeqIO.read(path, "genbank").seq).upper()


@pytest.fixture(scope="session")
def puc19() -> str:
    """pUC19, 2686 bp, circular plasmid (gx/corpus/inputs, see MANIFEST.json)."""
    text = (INPUTS / "puc19_L09137.2.txt").read_text().strip()
    assert len(text) == 2686 and set(text) <= set("ACGT")
    return text


@pytest.fixture(scope="session")
def jcv() -> str:
    """NC_001699.1, JC polyomavirus, 5130 bp, circular genome."""
    seq = _genbank_sequence(INPUTS / "NC_001699.1.gb")
    assert len(seq) == 5130
    return seq


@pytest.fixture(scope="session")
def bfdv() -> str:
    """AF071878.1, beak and feather disease virus, 1993 bp, circular genome."""
    seq = _genbank_sequence(INPUTS / "AF071878.1.gb")
    assert len(seq) == 1993
    return seq


@pytest.fixture(scope="session")
def col1a1() -> str:
    """NG_007400.1, human COL1A1 RefSeqGene, 24,544 bp, linear."""
    seq = _genbank_sequence(INPUTS / "NG_007400.1.gb")
    assert len(seq) == 24544 and set(seq) <= set("ACGT")
    return seq


@pytest.fixture(scope="session")
def mito_record():
    """NC_012920.1, the human mitochondrial rCRS: 16,569 bp, circular, table 2.

    The whole ``SeqRecord``, not just the bases -- the round needs its 13 CDS
    features, their ``/transl_table=2``, their ``/transl_except`` incomplete stop
    codons and above all their ``/translation`` qualifiers, which are recorded
    literals from NCBI and the oracle for the alternative-code requirements.
    """
    rec = SeqIO.read(INPUTS / "NC_012920.1.gb", "genbank")
    assert len(rec.seq) == 16569
    assert rec.annotations.get("topology") == "circular"
    return rec


@pytest.fixture(scope="session")
def mito_cds(mito_record):
    """The 13 rCRS CDS features as the APP reads them, keyed by gene name.

    Parsed with ``backend.formats.genbank.read_genbank`` rather than Biopython
    directly, because the subject of OR-2/OR-13/OR-14 is
    ``backend.services.annotations.restated_translation``, which takes the app's
    own ``Annotation`` model. Returns ``(sequence_str, {gene: Annotation})``.
    """
    from backend.formats.genbank import read_genbank

    with open(INPUTS / "NC_012920.1.gb") as handle:
        seq = read_genbank(handle)[0]
    assert seq.is_circular and len(seq.seq) == 16569
    by_gene = {
        str(a.qualifiers.get("gene")): a
        for a in seq.annotations
        if a.feature_type == "CDS"
    }
    assert len(by_gene) == 13, sorted(by_gene)
    return seq.seq, by_gene


@pytest.fixture(scope="session")
def corpus(puc19, jcv, bfdv, col1a1, mito_record) -> dict[str, str]:
    """Every pinned full-length nucleotide input, by the name the register uses."""
    return {
        "puc19_L09137.2": puc19,
        "NC_001699.1": jcv,
        "AF071878.1": bfdv,
        "NG_007400.1": col1a1,
        "NC_012920.1": str(mito_record.seq).upper(),
    }


@pytest.fixture
def find_orfs():
    """The subject's ORF finder, called the way an in-process caller calls it."""
    from backend.routers.analysis import find_orfs as _find_orfs

    return _find_orfs


@pytest.fixture
def translate_op():
    """``/sequences/{id}/edit`` with ``operation='translate'``, on a fresh record.

    Returns ``(result, warnings)``: the router's own ``SequenceEditResult`` and
    every warning raised while it ran, because whether the app handles or passes
    through Biopython's ``Partial codon`` warning is one of the things under test.
    """
    import warnings as _warnings

    from backend.models.sequence import MoleculeType, Sequence, SequenceEditRequest
    from backend.routers.sequences import _store, edit_sequence

    made: list[str] = []

    def _translate(seq: str, molecule_type: MoleculeType = MoleculeType.DNA):
        seq_id = f"gx-orf-{len(made)}"
        made.append(seq_id)
        _store[seq_id] = Sequence(
            id=seq_id, name=seq_id, seq=seq, molecule_type=molecule_type
        )
        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            result = edit_sequence(seq_id, SequenceEditRequest(operation="translate"))
        return result, list(caught)

    yield _translate

    from backend.routers.sequences import _store as store
    for seq_id in made:
        store.pop(seq_id, None)
        store.pop(f"{seq_id}_translate", None)
