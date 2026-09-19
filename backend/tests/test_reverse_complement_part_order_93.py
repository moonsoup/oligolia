"""#93: Reverse Complement reordered a spliced CDS's exons, changing the protein.

`flip_annotations` mapped every interval correctly — `[s, e)` -> `[L-e, L-s)`,
strand inverted — and then called `parts.reverse()`. A `CompoundLocation`'s parts
are already stored in the feature's own 5'-to-3' reading order, not in ascending
coordinate order, so reflecting each interval *already* preserves that order and
the extra reversal undid it. The record that came out was the right length with
every base still inside the feature; only the exon order was wrong, so the CDS
encoded a different protein and nothing signalled it.

The flip-flip round trip in `test_annotation_edits.py` cannot see this: reversing
a list twice restores it whether or not reversing it was right. So every test
here looks at **one** flip, and the oracle is
`Bio.SeqRecord.reverse_complement(features=True)` run on the same bytes plus the
translation of the CDS those bases actually spell.
"""

from __future__ import annotations

from io import StringIO

from Bio import SeqIO
from Bio.Seq import Seq

from backend.formats.genbank import read_genbank, write_genbank
from backend.services.annotations import flip_annotations

from ._spliced_circular_93 import (
    ORIGIN_SPANNING,
    SPLICED_CDS,
    SPLICED_CIRCULAR,
    SPLICED_PROTEIN,
)


def _app_record():
    return read_genbank(StringIO(SPLICED_CIRCULAR))[0]


def _bio_record():
    return SeqIO.read(StringIO(SPLICED_CIRCULAR), "genbank")


def _flipped(seq):
    """Reverse-complement a record exactly as `_commit_edit(flip=True)` does."""
    rc = str(Seq(seq.seq).reverse_complement())
    return seq.model_copy(update={
        "seq": rc,
        "annotations": flip_annotations(seq.annotations, len(seq.seq)),
        "length": len(rc),
    })


def _profile(record):
    """(type, parts, strand) per feature, from a Biopython record."""
    return [
        (f.type,
         [(int(p.start), int(p.end)) for p in f.location.parts],
         {1: "+", -1: "-"}.get(f.location.strand, "."))
        for f in record.features
    ]


def _app_profile(seq):
    return [
        (a.feature_type, [tuple(p) for p in a.parts], a.strand.value)
        for a in seq.annotations
    ]


def _feature(record, ftype: str):
    return next(f for f in record.features if f.type == ftype)


def _location_lines(genbank_text: str, key: str) -> list[str]:
    """Every location string written for feature key `key`."""
    return [
        line[21:].strip()
        for line in genbank_text.splitlines()
        if line.startswith(f"     {key:<16}")
    ]


def _saved(seq):
    """The record the app would write, parsed back with Biopython."""
    return SeqIO.read(StringIO(write_genbank([seq])), "genbank")


# ── the oracle: Biopython flips a record's features itself ────────────────────

def test_one_flip_matches_biopythons_own_reverse_complement() -> None:
    """A single flip, not a round trip — two reversals would cancel out.

    Sorted, so neither list order nor a lucky index can hide a difference:
    Biopython re-sorts the features it flips, the app keeps file order.
    """
    flipped = _flipped(_app_record())
    want = sorted(_profile(_bio_record().reverse_complement(features=True)))
    assert sorted(_app_profile(flipped)) == want


def test_every_feature_still_reads_the_same_bases_after_one_flip() -> None:
    """Flipping a molecule end-for-end does not edit a gene: the same bases are
    still in it, in the same order, read from the other strand."""
    before = _bio_record()
    after = _saved(_flipped(_app_record()))
    assert [str(f.extract(after.seq)) for f in after.features] == [
        str(f.extract(before.seq)) for f in before.features
    ]


# ── the consequence a user would see ─────────────────────────────────────────

def test_the_spliced_cds_still_starts_at_its_start_codon() -> None:
    saved = _saved(_flipped(_app_record()))
    bases = str(_feature(saved, "CDS").extract(saved.seq))
    assert bases.startswith("ATG"), bases[:12]
    assert bases == SPLICED_CDS
    assert str(Seq(bases).translate()) == SPLICED_PROTEIN


def test_the_origin_spanning_feature_keeps_its_bases_in_order() -> None:
    """The same bases in a different order is the whole defect, so assert the
    string and not just its length or its content as a set."""
    saved = _saved(_flipped(_app_record()))
    origin = _feature(saved, "rep_origin")
    assert str(origin.extract(saved.seq)) == ORIGIN_SPANNING


def test_the_exported_join_lists_its_exons_in_reading_order() -> None:
    """The error reaches the saved file, not just memory.

    Expected string comes from Biopython writing its *own* flipped record, so the
    INSDC convention is not restated by hand here — whatever `complement(join(…))`
    Biopython emits for the same molecule is what the app has to emit.
    """
    rc = _bio_record().reverse_complement(features=True)
    rc.annotations["molecule_type"] = "DNA"   # reverse_complement() does not carry it
    oracle = StringIO()
    SeqIO.write([rc], oracle, "genbank")
    assert _location_lines(write_genbank([_flipped(_app_record())]), "CDS") == \
        _location_lines(oracle.getvalue(), "CDS")


# ── nothing is lost ──────────────────────────────────────────────────────────

def test_the_flip_keeps_every_feature_and_the_length() -> None:
    original = _app_record()
    flipped = _flipped(original)
    assert len(flipped.annotations) == len(original.annotations) == 3
    assert flipped.length == original.length == 120
    assert flipped.seq == str(Seq(original.seq).reverse_complement())
