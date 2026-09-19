One root cause: the "contains the edit" branch in `shift_annotations` (`backend/services/annotations.py:74-76`) is reached by edits it was not written for, and it keeps the feature in both of them — collapsed to nothing in one case, still asserting a protein that no longer translates in the other.

Filed as one issue because both failures are the same branch taking the same wrong turn; the fix is a decision about what that branch should do.

Found by the Rockin-Robin `oligolia-gx-tests` contract, coordinates round 1, against commit `4bb77b3`.

## Requirements

**CO-10** — *A feature whose bases are entirely deleted must be reported as lost, not kept. Deleting exactly [5073, 5090) from NC_001699.1 — the whole of `repeat_region 5074..5090` — leaves no bases for that feature to describe, so it must appear in the dropped list. It must not survive as a zero-length annotation, which the exporter renders as the INSDC between-position `5073^5074`, a location that per the Feature Table 'points to a site between bases' and asserts something the original record never said.*

**CO-11** — *A feature that spans an edit is either adjusted so everything it asserts is still true, or reported as lost — never kept with a stale assertion. Deleting one base inside NC_001699.1's `CDS 526..1560` changes that CDS's own bases, so either the feature is dropped, or its `/translation` qualifier is updated; keeping `526..1559` (1034 bases, not a multiple of three) with the original 344-residue `/translation` exports a record whose stated protein does not translate from its stated bases.*

## Root cause

`backend/services/annotations.py:74-76`:

```python
elif p_start <= start and p_end >= end:
    # Contains the edit: keep the start, move the end.
    new_parts.append((p_start, p_end + delta))
```

This branch is tested **before** the overlap branch at `:77-80` that sets `lost = True`. Two kinds of edit reach it that should not:

1. **An exact-span delete.** `p_start <= start and p_end >= end` is satisfied exactly when the delete range equals the feature range, so the feature never reaches the overlap branch. The branch then keeps the start and moves the end, collapsing the interval to zero length.
2. **Any edit inside a feature.** The branch is right about coordinates — the function's docstring at `:52-54` says *"A feature that merely contains the edit grows or shrinks with it: those bases are still its bases"* — and wrong about everything else the feature asserts. `/translation` is a recorded claim about exactly those bases, and it is carried across untouched by `_rebuilt`.

The docstring at `:47-50` states the intended contract, which is not what happens:

> A feature whose bases were themselves edited is **dropped**, not truncated. Truncating would keep a feature whose sequence no longer matches what it claims to describe, which is the same class of quiet wrongness as the original bug — so the caller gets the list and can say which features went.

## Differential evidence

Input: `gx/corpus/inputs/NC_001699.1.gb` (JC polyomavirus, 5130 bp, circular), sha256 `6935cc5c717c934b1504a13d8a4c2f457c80adfbe7e9d0c16f1971a2d67a7dd0`.

### CO-10 — an exactly-deleted feature is kept as a zero-length between-position

| delete `[start, end)` | feature | `dropped` | kept parts | exported location |
|---|---|---|---|---|
| `[5073, 5090)` | `repeat_region 5074..5090` (#26) | *empty* | `[(5073, 5073)]` | `repeat_region   5073^5074` |
| `[492, 521)` | `intron 493..521` (#9) | *empty* | `[(492, 492)]` | `intron          492^493` |
| `[0, 12)` | `rep_origin 1..12` (#2) | target not dropped | `[(0, 0)]` | `rep_origin      0^1` |

Expected in every row: the feature appears in `shift_annotations`' `dropped` list and the user is told it went.

Why it matters beyond an in-memory oddity: INSDC's `n^n+1` means *"a site between bases n and n+1"* — a legitimate location type used for insertion points — so any reader will accept it without complaint and render a former 17-base repeat as a zero-width tick. The record now asserts something the original never said.

**`0^1` is additionally a malformed file**, independent of what one thinks should happen to a zero-length feature: a GenBank location counts from 1, so there is no base 0. Deleting NC_001699.1's first 12 bases exports both

```
rep_origin      0^1
rep_origin      join(5106..5118,0^1)
```

— so a circular plasmid trimmed at its origin exports a join with an invalid part. Biopython's own writer produced those strings without objection, so nothing downstream of the app will catch it either.

Reach: a user deleting a highlighted feature — the obvious way to remove one — hits this on the first try.

### CO-11 — an edited CDS keeps its old `/translation`

Edit: delete one base at `[700, 701)`, inside `CDS 526..1560` (feature #11), `delta = -1`.

| | expected | actual |
|---|---|---|
| CDS fate | in `dropped`, **or** kept with `/translation` equal to `str(extract.translate(table=1, to_stop=True))` of the edited bases | kept as `526..1559`, absent from `dropped` |
| length | a whole number of codons | **1034** bases (1035 before), not a multiple of three; Biopython emits `Partial codon, len(sequence) not a multiple of three` when asked to translate it |
| `/translation` | restated | unchanged 344-residue `MGAALALLGDLVATVSEAAAATGFSVAEIA…` — the bases now encode a different protein from the frameshift onward |

The frame check is a convention-free second angle: whatever the right policy is for a frameshifted CDS, a complete (non-partial, not originally short) CDS that is not a multiple of three is not a valid annotation, and the only thing that made it invalid is the edit.

**Equal-length replace — the clearest case, no coordinate arithmetic to argue about.** Replacing `[1700, 1710)` inside `CDS 1469..2533` (feature #15) with ten bases, `delta = 0`:

- The CDS is kept at `1469..2533`, every coordinate unchanged, `/translation` unchanged.
- The bases now encode `…FESDSPNRDMLPCYSVARIPLPN…`; the record still claims `…FESDSPNRDMQKKNSVARIPLPN…`.

Nothing moved, and the record is still wrong. An equal-length replacement inside a feature is the ordinary way a user applies a point mutation or a codon swap — the app will export the edited plasmid still advertising the unmutated protein.

The GUI reports both edits as clean: the "… annotation(s) removed" message at `gui/panels/sequence_panel.py:951-959` is driven by the `dropped` list, which is empty.

## Relationship to closed issues

This is the **direct continuation of #59** ("Sequence edits never shift annotations"). #59's fix gave the GUI `shift_annotations`, and the shifting itself is sound — CO-8 (an edit elsewhere leaves an untouched feature's bases byte-identical, for all four pinned edits) and CO-9 (downstream parts shift by exactly `delta`, upstream parts do not move, re-derived by index from the raw sequence) both **hold**. What #59 left unfixed is the case its own docstring names: features now move, but a feature whose own bases were edited is still kept and still says something false.

## Pinned source anchors

| anchor | sha256 |
|---|---|
| `insdc/FT_current.txt#123^124                   Points to a site between bases 123 and 124` | `14ea7e30324fa7bcb50b108bce45cf80de2eb37075d06be3f36ffaf9f7245444` |
| `insdc/FT_current.txt#join(12..78,134..202)     Regions 12 to 78 and 134 to 202 should be joined to form` | `14ea7e30324fa7bcb50b108bce45cf80de2eb37075d06be3f36ffaf9f7245444` |
| `biopython-1.85/Bio/SeqFeature.py#    def extract(self, parent_sequence, references=None):` | `9bdb8d0386de1848092e8c7901e80b0955074ff0cafc7e89c4e41ae4d07b7d56` |
| `inputs/NC_001699.1.gb#     repeat_region   5074..5090` | `6935cc5c717c934b1504a13d8a4c2f457c80adfbe7e9d0c16f1971a2d67a7dd0` |
| `inputs/NC_001699.1.gb#     CDS             526..1560` | `6935cc5c717c934b1504a13d8a4c2f457c80adfbe7e9d0c16f1971a2d67a7dd0` |

## Sealed failing tests

Seal digest: `2ef985275ccf556d2f29023f030d029fabb0a5e5d723ad3c33591a33765987f9` (`gx/coords/seal.json`). Same 42 failed / 158 passed on Linux and macOS. **The sealed tests must not be edited.**

```
gx/coords/test_coords_edit_shift.py::test_co_10_a_wholly_deleted_feature_is_reported_not_kept_empty[exact_delete_of_repeat_region_5074..5090]
gx/coords/test_coords_edit_shift.py::test_co_10_a_wholly_deleted_feature_is_reported_not_kept_empty[exact_delete_of_intron_493..521]
gx/coords/test_coords_edit_shift.py::test_co_10_a_wholly_deleted_feature_is_reported_not_kept_empty[exact_delete_of_rep_origin_1..12]
gx/coords/test_coords_edit_shift.py::test_co_10_no_zero_length_between_position_is_exported[exact_delete_of_repeat_region_5074..5090]
gx/coords/test_coords_edit_shift.py::test_co_10_no_zero_length_between_position_is_exported[exact_delete_of_intron_493..521]
gx/coords/test_coords_edit_shift.py::test_co_10_no_zero_length_between_position_is_exported[exact_delete_of_rep_origin_1..12]
gx/coords/test_coords_edit_shift.py::test_co_10_no_exported_location_names_base_zero
gx/coords/test_coords_edit_shift.py::test_co_11_a_cds_whose_bases_changed_does_not_keep_a_stale_translation
gx/coords/test_coords_edit_shift.py::test_co_11_the_shortened_cds_is_not_left_out_of_frame
gx/coords/test_coords_edit_shift.py::test_co_11_a_same_length_replace_inside_a_cds_is_not_silent
```

## Acceptance

- The ten sealed gx tests listed above pass, unedited.
- `make check` stays green — including the existing `shift_annotations` coverage in `backend/tests/test_annotation_edits.py`.
- Nothing is removed. In particular CO-8 and CO-9 keep holding: an edit elsewhere must still leave an untouched feature's bases byte-identical, and a feature strictly after the edit must still shift by exactly `delta` while one strictly before it does not move. Whichever policy is chosen for CO-11 (drop the feature, or restate `/translation`), features that merely *contain* an edit without their own bases changing — which is impossible for an interior edit but applies to a delete flush against a boundary — must not start disappearing.
