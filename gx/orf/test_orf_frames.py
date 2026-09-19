"""OR-3, OR-4, OR-5, OR-7, OR-8 — six frames, both strands, by index.

Every assertion here is over full-length pinned records, and every expected
value comes either from ``Bio.Seq.translate`` on the same bytes, from gc.prt, or
from re-extracting the subject's own reported interval out of the input. The
subject's output is never its own oracle.
"""

from __future__ import annotations

import pytest
from orf_oracle import (
    biopython_translation,
    find_orfs_oracle,
    frame_substring,
    ncbi_amino_acids,
    orf_bases,
    top_strand_span,
)

INPUTS = ["puc19_L09137.2", "NC_001699.1", "AF071878.1", "NG_007400.1", "NC_012920.1"]
STOPS = {c for c, a in ncbi_amino_acids(1).items() if a == "*"}


@pytest.mark.parametrize("name", INPUTS)
def test_or_3_each_frame_agrees_with_biopython_translate(name, corpus, find_orfs) -> None:
    """An ORF in frame f is the segment of Bio.Seq.translate(frame f) it claims.

    Not 'a protein of the right length': the ORF's coordinates fix a residue
    offset into the frame's translation, and the reported protein must be
    exactly the slice at that offset.
    """
    sequence = corpus[name]
    n = len(sequence)
    frame_protein = {f: biopython_translation(frame_substring(sequence, f), 1) for f in (1, 2, 3, -1, -2, -3)}

    for orf in find_orfs(sequence, min_length_aa=30).orfs:
        # Where the ORF starts within its own frame's strand, in residues.
        s, e = top_strand_span(n, orf.start, orf.end, orf.frame)
        offset = abs(orf.frame) - 1
        assert (s - offset) % 3 == 0, f"{name} {orf.frame} [{orf.start}:{orf.end}) is not in frame"
        first = (s - offset) // 3
        span_aa = (e - s) // 3
        segment = frame_protein[orf.frame][first : first + span_aa]
        # The subject renders the initiator M by convention; compare past it.
        assert segment[1:].rstrip("*") == orf.protein[1:], (
            f"{name} frame {orf.frame} [{orf.start}:{orf.end}): reported protein "
            f"is not Bio.Seq.translate's residues {first}..{first + span_aa} of that frame"
        )


@pytest.mark.parametrize("name", INPUTS)
def test_or_4_reverse_strand_coordinates_map_back_to_the_top_strand(name, corpus, find_orfs) -> None:
    """A minus-frame ORF re-extracts from the TOP strand at [start, end)."""
    sequence = corpus[name]
    n = len(sequence)
    minus = [o for o in find_orfs(sequence, min_length_aa=30).orfs if o.frame < 0]
    assert minus, f"{name} reports no minus-strand ORF to check"

    for orf in minus:
        bases = orf_bases(sequence, orf.start, orf.end, orf.frame)
        assert bases[:3] == orf.start_codon, (
            f"{name} frame {orf.frame} [{orf.start}:{orf.end}): reverse complement "
            f"of the reported span starts {bases[:3]}, reported start_codon "
            f"{orf.start_codon}"
        )
        assert biopython_translation(bases, 1)[1:].rstrip("*") == orf.protein[1:]
        # The mapping is its own inverse, so the span must round-trip.
        s, e = top_strand_span(n, orf.start, orf.end, orf.frame)
        assert top_strand_span(n, s, e, orf.frame) == (orf.start, orf.end)


@pytest.mark.parametrize("name", INPUTS)
def test_or_5_every_orf_re_extracts_start_to_stop(name, corpus, find_orfs) -> None:
    """Start codon at the 5' end, stop codon at the 3' end, none in between."""
    sequence = corpus[name]
    for orf in find_orfs(sequence, min_length_aa=30).orfs:
        bases = orf_bases(sequence, orf.start, orf.end, orf.frame)
        where = f"{name} frame {orf.frame} [{orf.start}:{orf.end})"
        assert bases[:3] == "ATG", f"{where}: opens on {bases[:3]!r}, not a default start codon"
        assert bases[:3] == orf.start_codon, f"{where}: start_codon field disagrees with the bases"
        codons = [bases[i : i + 3] for i in range(0, len(bases) // 3 * 3, 3)]
        interior = [(i, c) for i, c in enumerate(codons[:-1]) if c in STOPS]
        assert interior == [], f"{where}: stop codon inside the ORF at codon(s) {interior[:5]}"
        assert codons[-1] in STOPS, (
            f"{where}: ends on {codons[-1]!r}, which is not a table-1 stop codon, "
            f"and the result carries no field saying the ORF is partial"
        )


@pytest.mark.parametrize("name", INPUTS)
@pytest.mark.parametrize("min_length_aa", [30, 100])
def test_or_7_min_length_is_honoured_in_both_directions(
    name, min_length_aa, corpus, find_orfs
) -> None:
    """Set equality against a second ORF finder written from gc.prt."""
    sequence = corpus[name]
    result = find_orfs(sequence, min_length_aa=min_length_aa)

    short = [(o.frame, o.start, o.end, o.length_aa) for o in result.orfs if o.length_aa < min_length_aa]
    assert short == [], f"{name}: ORFs below min_length_aa={min_length_aa}: {short[:5]}"

    reported = {(o.frame, o.start, o.end, o.length_aa, o.protein) for o in result.orfs}
    expected = {
        (o["frame"], o["start"], o["end"], o["length_aa"], "M" + o["protein"][1:])
        for o in find_orfs_oracle(
            sequence, min_length_aa=min_length_aa, table_id=1, start_codons={"ATG"}
        )
    }
    assert reported == expected, (
        f"{name} m={min_length_aa}: "
        f"{len(expected - reported)} ORF(s) missed, {len(reported - expected)} spurious"
    )
    assert result.sequence_length == len(sequence)
    assert result.total_found == len(result.orfs)


@pytest.mark.parametrize("name", INPUTS)
def test_or_8_every_orf_is_a_whole_number_of_codons(name, corpus, find_orfs) -> None:
    """``end - start == length_nt`` and ``length_nt % 3 == 0``."""
    for orf in find_orfs(corpus[name], min_length_aa=30).orfs:
        where = f"{name} frame {orf.frame} [{orf.start}:{orf.end})"
        assert orf.end - orf.start == orf.length_nt, f"{where}: length_nt={orf.length_nt}"
        assert orf.length_nt % 3 == 0, (
            f"{where}: length_nt={orf.length_nt} is not a multiple of three — the "
            f"reported interval ends part-way through a codon"
        )
