**Severity: high — silent data corruption in a routine one-click operation.** The record that comes out is well-formed, the same length, with every base still inside the feature and the correct strand; only the order of a spliced feature's exons is wrong, so the CDS encodes a different protein and nothing in the app or in the file signals it.

Found by the Rockin-Robin `oligolia-gx-tests` contract, coordinates round 1, against commit `4bb77b3`.

## Requirement

**CO-12** — *Reverse-complementing a record maps every feature onto the new strand without changing what the feature reads: for every feature, `extract` on the reverse-complemented record must return exactly the string `extract` returned on the original record, because reverse-complementing a molecule does not change which bases a gene is made of or the order it reads them in. This must hold for compound features too — a spliced CDS's exons must still be listed in 5'-to-3' feature order, which for a minus-strand CompoundLocation means descending coordinates.*

## Root cause

`backend/services/annotations.py:104-105`

```python
parts = [(seq_len - p_end, seq_len - p_start) for p_start, p_end in _parts_of(ann)]
parts.reverse()
```

The interval map on line 104 is correct. The `parts.reverse()` on line 105 is not.

A `CompoundLocation`'s parts are already stored in the feature's own 5'-to-3' reading order, not in ascending coordinate order. Biopython's INSDC writer states the convention outright:

> we expect the CompoundLocation and its parts to all be marked as strand == -1, and to be in the order 19:100 then 0:10.
> — `gx/corpus/biopython-1.85/Bio/SeqIO/InsdcIO.py:334-339`

Reflecting each interval through `L - x` already preserves that reading order. Reversing the list afterwards undoes it.

Every single-part feature is unaffected, which is why 26 of NC_001699.1's 28 features pass and why this survived: it only shows on a record that has a `join()`, and only under Reverse Complement.

## Why it has not been caught

`flip_annotations`' own docstring (`backend/services/annotations.py:98-99`) says:

> Applying this twice returns the original, which is how the mapping is checked rather than reasoned about.

**A flip-flip round trip cannot detect this defect.** Reversing a list twice restores it whether or not reversing it was correct in the first place, so `backend/tests/test_annotation_edits.py:152` (`test_flipping_twice_returns_the_original`) passes with the bug present and would pass with it fixed. The extracted bases, not the round trip, are what pins this.

The neighbouring `backend/tests/test_annotation_edits.py:144-149` (`test_reverse_complement_reverses_the_part_order`) actively locks the wrong behaviour in: it flips a synthetic plus-strand feature with parts `[(4, 10), (29, 40)]` and asserts the result is `[(60, 71), (90, 96)]` — ascending, on what is now a minus-strand feature. By the convention above the correct result is `[(90, 96), (60, 71)]`. **Fixing this requires updating that expectation**; it is the only place in `make check` that asserts the current behaviour.

## Differential evidence

Input: `gx/corpus/inputs/NC_001699.1.gb` (JC polyomavirus, 5130 bp, circular), sha256 `6935cc5c717c934b1504a13d8a4c2f457c80adfbe7e9d0c16f1971a2d67a7dd0`.

Oracle: `Bio.SeqRecord.reverse_complement(features=True)` from Biopython 1.85, run on the same bytes. It agrees with the app on every interval and every strand, and disagrees only on the order of the two compound part lists.

| feature | oracle (Biopython) | app (`flip_annotations`) |
|---|---|---|
| `CDS complement(join(2603..4426,4771..5013))`, index 20 | `[(117, 360), (704, 2528)]` strand `+` | `[(704, 2528), (117, 360)]` strand `+` |
| `rep_origin join(5118..5130,1..12)`, index 1 | `[(0, 13), (5118, 5130)]` strand `-` | `[(5118, 5130), (0, 13)]` strand `-` |

Exported locations:

| | expected | actual |
|---|---|---|
| CDS | `join(118..360,705..2528)` | `join(705..2528,118..360)` |
| rep_origin | `complement(join(5119..5130,1..13))` | `complement(join(1..13,5119..5130))` |

`join(705..2528,118..360)` is a descending join on the plus strand, which is only meaningful for a feature that crosses the origin. This one does not cross it, so the location is both wrong and misleading about why.

Extracted bases — the consequence:

- Large T antigen CDS, 2067 bases either way.
  - expected: `ATGGACAAAGTGCTGAATAGGGAGGAATCCATGGAGCTTATGGATTTATTAGGCCTTGAT`
  - actual: `GTGCCAACCTATGGAACAGATGAATGGGAATCCTGGTGGAATACATTTAATGAGAAGTGG` — the second exon read first, so **the CDS no longer begins with a start codon**.
- Origin-spanning `rep_origin`, 25 bases either way.
  - expected: `GGAGGCGGAGGCGGCCTCGGCCTCC`
  - actual: `GCCTCGGCCTCCGGAGGCGGAGGCG` — the same 25 bases with the two parts swapped.

2 of 28 features differ; both are the compound ones. The error reaches the saved file, not just memory.

Driven end-to-end through the shipped widget (`gui/panels/sequence_panel.py:_apply_op` → `_commit_edit(flip=True)` at `:958-961`, saved via `gui/main_window.py:210` → `write_genbank`), the same two features come back wrong: `join{[5118:5130](-), [0:13](-)}` and `join{[704:2528](+), [117:360](+)}`. The companion assertion that the feature table still lists all 28 rows **passes** — the app gives the user no signal that anything changed.

## Relationship to closed issues

Not a regression of #57 or #59, and not an unfixed half of either — `flip_annotations` was added for #59 ("Annotations are never shifted … Insert, delete, replace and reverse-complement all leave feature coordinates untouched") and the coordinates it produces are right. The part *ordering* is a defect introduced with that fix, and #59's own tests could not see it for the flip-flip reason above.

## Pinned source anchors

| anchor | sha256 |
|---|---|
| `biopython-1.85/Bio/SeqIO/InsdcIO.py#strand == -1, and to be in the order 19:100 then 0:10.` | `2c12bd14c507928e931b17bcd5ece6a44175479aa17e9efd7b1fec40d8b770a9` |
| `biopython-1.85/Bio/SeqFeature.py#    def extract(self, parent_sequence, references=None):` | `9bdb8d0386de1848092e8c7901e80b0955074ff0cafc7e89c4e41ae4d07b7d56` |
| `insdc/FT_current.txt#complement(location)` | `14ea7e30324fa7bcb50b108bce45cf80de2eb37075d06be3f36ffaf9f7245444` |
| `inputs/NC_001699.1.gb#     CDS             complement(join(2603..4426,4771..5013))` | `6935cc5c717c934b1504a13d8a4c2f457c80adfbe7e9d0c16f1971a2d67a7dd0` |

## Sealed failing tests

Seal digest: `2ef985275ccf556d2f29023f030d029fabb0a5e5d723ad3c33591a33765987f9` (`gx/coords/seal.json`). Same 42 failed / 158 passed on Linux and macOS. **The sealed tests must not be edited.**

```
gx/coords/test_coords_reverse_complement.py::test_co_12_every_feature_still_extracts_the_same_string
gx/coords/test_coords_reverse_complement.py::test_co_12_matches_biopythons_own_reverse_complement
gx/coords/test_coords_reverse_complement.py::test_co_12_compound_parts_stay_in_reading_order[20-CDS]
gx/coords/test_coords_reverse_complement.py::test_co_12_compound_parts_stay_in_reading_order[1-rep_origin]
gx/coords/test_coords_reverse_complement.py::test_co_12_compound_exports_with_its_exons_in_the_right_order[20-CDS]
gx/coords/test_coords_reverse_complement.py::test_co_12_compound_exports_with_its_exons_in_the_right_order[1-rep_origin]
gx/coords/test_coords_reverse_complement.py::test_co_12_spliced_cds_still_starts_at_its_start_codon
gx/coords/test_coords_reverse_complement.py::test_co_12_origin_spanning_feature_keeps_its_bases
gx/coords/test_coords_reverse_complement.py::test_co_12_the_same_holds_after_export_and_re_read
gx/coords/test_coords_gui_edit.py::test_co_12_reverse_complement_through_the_panel_keeps_every_feature_s_bases
```

## Acceptance

- The ten sealed gx tests listed above pass, unedited.
- `make check` stays green. The expectation in `backend/tests/test_annotation_edits.py:144-149` must be corrected rather than deleted — it currently asserts the wrong order.
- Nothing is removed: all 28 features of NC_001699.1 still survive the flip, `flip_annotations` still returns its input length, and the 158 currently-passing gx coords tests stay passing.
