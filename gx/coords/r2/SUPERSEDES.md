# Coordinates round 2 — supersedes `gx-co-53`

**Adjudication: TEST-DEFECT.** The app at `d22ffc8` ("Carry the annotations and the
topology through an API edit (#96)") is right. The sealed round-1 test is wrong.

## What this replaces

| | |
|---|---|
| Sealed test id | `gx-co-53` (`gx/coords/tests.json`) |
| Node id | `gx/coords/test_coords_edit_endpoint.py::test_co_14_surviving_annotations_point_at_the_right_bases` |
| Covers | CO-14 (and CO-9) |
| Replaced by | `gx-co2-1` … `gx-co2-7`, covering CO2-1, CO2-2, CO2-3 |

The sealed test stays on disk and keeps failing. That failure is the record of
what round 1 got wrong; it is not something to hide or edit away. Nothing sealed
has been touched — round 2 is additive, under `gx/coords/r2/`.

## The defect

The assertion is:

```python
before = {(f.type, str(f.location)): str(f.extract(jcv_bio.seq)) for f in jcv_bio.features}
...
wrong = [a.feature_type for a in stored.annotations
         if spliced(stored.seq, [tuple(p) for p in a.parts], a.strand.value)
         not in before.values()]
assert wrong == [], wrong
```

It requires **every** surviving annotation to extract a byte string the
**pre-edit** record already held. It carries no exemption for a feature whose
own bases the edit changed — no `touched()` guard, which is exactly the guard
its four sibling tests use for this purpose:

- `gx/coords/test_coords_edit_shift.py:60` (CO-8)
- `gx/coords/test_coords_edit_shift.py:76` (CO-8)
- `gx/coords/test_coords_gui_edit.py:84` (CO-8)
- `gx/coords/test_coords_gui_edit.py:105`, `:128` (CO-9)

The test's own docstring says "every feature **downstream** of the insert still
has to read its own bases" — downstream, which is CO2-2 here. The assertion it
actually makes is over all features, including the containing ones.

## The evidence

The pinned edit is CO-14's first: insert 30 `A`s at position 3000 of
NC_001699.1 (5130 bp, circular). Five of the record's 28 features have an
interval that strictly contains position 3000, and the app grows each of them by
exactly the inserted bases:

| # | feature | strand | parts before | parts after | extract |
|---|---|---|---|---|---|
| 0 | `source` | + | `[(0, 5130)]` | `[(0, 5160)]` | 5130 → 5160 bp |
| 16 | `prim_transcript` | − | `[(2526, 5115)]` | `[(2526, 5145)]` | 2589 → 2619 bp |
| 19 | `gene` | − | `[(2602, 5013)]` | `[(2602, 5043)]` | 2411 → 2441 bp |
| 20 | `CDS` | − | `[(4770, 5013), (2602, 4426)]` | `[(4800, 5043), (2602, 4456)]` | 2067 → 2097 bp |
| 21 | `exon` | − | `[(2602, 4426)]` | `[(2602, 4456)]` | 1824 → 1854 bp |

In each case the post-edit extract is the feature's own pre-edit extract with
the 30 inserted bases spliced in at the corresponding offset — reverse
complemented for the four minus-strand features, and for feature 20 landing
inside the second part of the join, the one that contains position 3000. Those
are precisely the five names the sealed test reports:

```
AssertionError: ['source', 'prim_transcript', 'gene', 'CDS', 'exon']
gx/coords/test_coords_edit_endpoint.py:95
1 failed, 13 passed
```

Nothing else about CO-14 is in question: the other 13 assertions in the same
file pass, and the other 23 features extract byte-identically.

**1. The assertion contradicts CO-14's own statement.** CO-14 requires the
stored record to carry the record's features "shifted for the edit, **the way
the GUI's own edit path does via `shift_annotations`**". That function's
contract, in its own docstring at `backend/services/annotations.py`, is:

> A feature that merely *contains* the edit grows or shrinks with it: those
> bases are still its bases. A feature strictly after the edit shifts by the net
> length change; one strictly before it does not move.

A feature that grew by 30 bases cannot extract a byte string the pre-edit record
held. The only way to satisfy the sealed assertion is to drop all five — which
CO-14's own wording forbids.

**2. Passing it would reintroduce the bug #96 fixes.** `source` is
`1..5130`, the whole plasmid. Dropping it on every insertion is the silent
annotation loss #96 exists to fix.

**3. Passing it would turn sealed round-1 tests red.** Simulated by dropping the
five containing features and re-running the sealed GUI oracle's own logic: the
saved record comes back with 23 features instead of 28, so
`test_co_9_the_saved_file_after_an_insert_still_describes_the_same_features`
(`gx/coords/test_coords_gui_edit.py:92`), which pairs `saved.features[i]` with
`before[i]` positionally, misaligns and reports 14 wrong features, beginning:

```
(1, 'rep_origin', '[0:12](+)')
(2, 'repeat_region', '[10:109](+)')
(3, 'repeat_region', '[109:207](+)')
```

## What round 2 asserts instead

CO-14, restated as the three clauses it actually contains
(`gx/coords/r2/register.json`):

- **CO2-1** — a feature that strictly contains an insert grows by the inserted
  bases and still describes its own feature: the containing part keeps its start
  and moves its end by `+30`, and its extract is its own pre-edit bases with the
  insert spliced in. All 28 features are stored; the topology survives.
- **CO2-2** — a part strictly downstream shifts by exactly `delta` and reads
  byte-identical bases; a feature no part of which contains the edit extracts
  byte-identically. Worked example: the origin-spanning
  `join(5118..5130,1..12)`, where one part moves and the other does not.
- **CO2-3** — a feature strictly upstream is untouched: same parts, same strand,
  same bases.

All three are checked across all three pinned splice edits where they apply
(insert 30 @ 3000; delete [2100, 2200); replace [1600, 1650) with 10 bases). No
feature of NC_001699.1 partially overlaps any of them, so no feature has grounds
to be dropped — the suite asserts that as a checked precondition rather than
assuming it.

## Discrimination

The corrected tests are not vacuous:

- at `d22ffc8` (#96 fixed): **21 passed**
- at `d22ffc8^` (#96 unfixed): **21 failed**

Expected values come from `Bio.SeqIO`/`SeqFeature.extract` on the pinned bytes,
from re-derivation by index on the raw string, and from literals in
`gx/coords/r2/register.json` — never from the app. Shared helpers live in
`coords_r2_oracle.py` and are imported by that name, so round 2 does not collide
with round 1's two `conftest.py` modules.
