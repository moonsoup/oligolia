Three losses, one root cause: `backend/formats/genbank.py` flattens a Biopython location to `(int, int)` pairs on the way in (`:32-37`) and rebuilds a plain `FeatureLocation` on the way out (`:95-96`). Everything a location string carries *besides* its two integers — the `<`/`>` position class, the `order` vs `join` operator, the remote accession on a part — is gone before the value is ever stored, and cannot be recovered at export.

Filed as one issue because all three are the same flattening in the same reader/writer pair; fixing it for one boundary fixes it for all three.

Found by the Rockin-Robin `oligolia-gx-tests` contract, coordinates round 1, against commit `4bb77b3`.

## Requirements

**CO-5** — *Partial boundaries survive the round trip. INSDC gives `<` and `>` a meaning no exact coordinate carries — `<1..888` means the feature starts before the first sequenced base — so a record whose file says `<108..1007` must not be exported as `108..1007`. The three partial boundaries in the pinned real records must still be partial in the subject's export.*

**CO-6** — *`order()` must not be exported as `join()`. INSDC defines them as different assertions — join means the elements 'should be joined (placed end-to-end) to form one contiguous sequence', order means only that they 'can be found in the specified order ... but nothing is implied about the reasonableness about joining them' — so rewriting one as the other changes what the record claims.*

**CO-7** — *A remote reference part — INSDC's `J00194.1:100..202`, 'bases 100 to 202, inclusive, in the entry ... with primary accession number J00194' — must be preserved or the feature rejected. It must never be exported as a bare local interval, because that silently reassigns the feature to this record's own bases 100..202.*

## Root cause

`backend/formats/genbank.py:32-37` (read):

```python
start = int(feat.location.start)
end = int(feat.location.end)
parts = [(int(p.start), int(p.end)) for p in feat.location.parts]
```

`backend/formats/genbank.py:95-96` (write):

```python
locs = [FeatureLocation(a, b, strand=strand_val) for a, b in parts]
location = locs[0] if len(locs) == 1 else CompoundLocation(locs)
```

- `int(p.start)` discards the position *class*. `BeforePosition` and `AfterPosition` are the only carriers of partiality, and both subclass `int`.
- `location.operator` is read nowhere, and `CompoundLocation`'s operator defaults to `"join"`.
- `p.ref` (the remote accession) is read nowhere, though Biopython does populate it on the parsed location — the information is available at the point it is thrown away.

## Differential evidence

### CO-5 — partial boundaries become exact

Inputs: `gx/corpus/inputs/AF071878.1.gb` (sha256 `f03eccea45a1b735d94141f8444af6990a1f6d4c109ae647f6f72442e0aa0e0f`) and `gx/corpus/inputs/NC_001699.1.gb` (sha256 `6935cc5c717c934b1504a13d8a4c2f457c80adfbe7e9d0c16f1971a2d67a7dd0`). Oracle: `Bio.SeqIO.read` on the pinned bytes vs on the app's own export.

| record / feature | file says | app exports |
|---|---|---|
| AF071878.1 #2 `gene` | `<108..1007` | `108..1007` |
| AF071878.1 #3 `CDS` | `<108..1007` | `108..1007` |
| AF071878.1 #6 `gene` | `complement(1240..>1974)` | `complement(1240..1974)` |
| AF071878.1 #7 `CDS` | `complement(1240..>1974)` | `complement(1240..1974)` |
| NC_001699.1 #21 `exon` | `complement(<2603..4426)` | `complement(2603..4426)` |

Re-read as Biopython locations: `[<107:1007](+)` in → `[107:1007](+)` out; `[<2602:4426](-)` in → `[2602:4426](-)` out. The literal `<` / `>` character is absent from the saved bytes.

Why this matters: AF071878.1's CDS is annotated `<108..1007` precisely because the coding sequence runs off the sequenced fragment. Re-saving it as `108..1007` asserts that base 108 is the A of the start codon — a claim the original explicitly declined to make. Downstream tools act on the difference: a `<` CDS is exempt from the start-codon and multiple-of-three checks an exact one must pass. Five of the eight non-source features in AF071878.1 are affected, i.e. this fires on the first real partial record a user opens.

**Any test comparing coordinates is blind to this**: `BeforePosition(5) == ExactPosition(5)` is `True` (`gx/corpus/biopython-1.85/Bio/SeqFeature.py:2209ff`), which is likely why it survived #57's round-trip work.

### CO-6 — `order()` exports as `join()`

Probe record: a minimal GenBank record over the first 300 bases of AF071878.1, so the sequence under test is real.

| in | out |
|---|---|
| `order(11..20,51..60)` | `join(11..20,51..60)` |
| `join(101..110,141..150)` (control) | `join(101..110,141..150)` — correct |

Parsed in as `order{[10:20](+), [50:60](+)}`. The control passing in the same record isolates the failure to the dropped operator rather than to compound handling in general — and means a fix must not be made by breaking `join`.

Lower reach than CO-5: `order()` is rare and no record in the pinned corpus uses it. But the record comes back asserting the parts splice together, which is the one thing `order` exists to avoid saying.

### CO-7 — a remote reference becomes local bases

Same probe record.

| in | out |
|---|---|
| `join(201..210,J00194.1:100..202)` | `join(201..210,100..202)` |

Parsed in as `join{[200:210](+), J00194.1[99:202](+)}` — the accession is right there on the part.

This is the worst-behaved of the three because it does not merely drop information, it **reassigns** it: bases 100..202 of *this* record are now claimed as part of the feature, and nothing in the output marks the substitution. Remote references are uncommon in the records this app targets, so severity is low; the failure mode is silent corruption rather than a visible gap.

## Relationship to closed issues

This is the **unfixed remainder of #57** ("GenBank round-trip collapses join() locations and stringifies multi-value qualifiers"). #57's fix direction was *"keep `CompoundLocation` parts through load, edit and export"*, and the parts themselves now do survive — the comments at `backend/formats/genbank.py:34-36` and `:91-93` both cite #57, and requirement CO-2 (the full coordinate round trip is the identity on both real records) and CO-3 (join parts survive, including the origin-spanning ones) both **hold**. What was not kept is everything about a location other than its integers. Same file, same two code paths, same class of loss.

## Pinned source anchors

| anchor | sha256 |
|---|---|
| `insdc/FT_current.txt#<1..888                   The feature starts before the first sequenced base and` | `14ea7e30324fa7bcb50b108bce45cf80de2eb37075d06be3f36ffaf9f7245444` |
| `insdc/FT_current.txt#1..>888                   The feature starts at the first sequenced base and` | `14ea7e30324fa7bcb50b108bce45cf80de2eb37075d06be3f36ffaf9f7245444` |
| `insdc/FT_current.txt#order(location,location, ... location)` | `14ea7e30324fa7bcb50b108bce45cf80de2eb37075d06be3f36ffaf9f7245444` |
| `insdc/FT_current.txt#join(1..100,J00194.1:100..202)` | `14ea7e30324fa7bcb50b108bce45cf80de2eb37075d06be3f36ffaf9f7245444` |
| `insdc/FT_current.txt#J00194.1:100..202         Points to bases 100 to 202, inclusive, in the entry (in` | `14ea7e30324fa7bcb50b108bce45cf80de2eb37075d06be3f36ffaf9f7245444` |
| `biopython-1.85/Bio/SeqFeature.py#class BeforePosition(int, Position):` | `9bdb8d0386de1848092e8c7901e80b0955074ff0cafc7e89c4e41ae4d07b7d56` |
| `biopython-1.85/Bio/SeqFeature.py#class AfterPosition(int, Position):` | `9bdb8d0386de1848092e8c7901e80b0955074ff0cafc7e89c4e41ae4d07b7d56` |
| `biopython-1.85/Bio/SeqFeature.py#class CompoundLocation(Location):` | `9bdb8d0386de1848092e8c7901e80b0955074ff0cafc7e89c4e41ae4d07b7d56` |
| `inputs/AF071878.1.gb#     CDS             complement(1240..>1974)` | `f03eccea45a1b735d94141f8444af6990a1f6d4c109ae647f6f72442e0aa0e0f` |
| `inputs/NC_001699.1.gb#     exon            complement(<2603..4426)` | `6935cc5c717c934b1504a13d8a4c2f457c80adfbe7e9d0c16f1971a2d67a7dd0` |

## Sealed failing tests

Seal digest: `2ef985275ccf556d2f29023f030d029fabb0a5e5d723ad3c33591a33765987f9` (`gx/coords/seal.json`). Same 42 failed / 158 passed on Linux and macOS. **The sealed tests must not be edited.**

```
gx/coords/test_coords_export.py::test_co_5_partial_boundaries_survive_export[AF071878.1.gb-2-gene]
gx/coords/test_coords_export.py::test_co_5_partial_boundaries_survive_export[AF071878.1.gb-3-CDS]
gx/coords/test_coords_export.py::test_co_5_partial_boundaries_survive_export[AF071878.1.gb-6-gene]
gx/coords/test_coords_export.py::test_co_5_partial_boundaries_survive_export[AF071878.1.gb-7-CDS]
gx/coords/test_coords_export.py::test_co_5_partial_boundaries_survive_export[NC_001699.1.gb-21-exon]
gx/coords/test_coords_export.py::test_co_5_exported_location_string_still_carries_the_angle_bracket[AF071878.1.gb-2-gene]
gx/coords/test_coords_export.py::test_co_5_exported_location_string_still_carries_the_angle_bracket[AF071878.1.gb-3-CDS]
gx/coords/test_coords_export.py::test_co_5_exported_location_string_still_carries_the_angle_bracket[AF071878.1.gb-6-gene]
gx/coords/test_coords_export.py::test_co_5_exported_location_string_still_carries_the_angle_bracket[AF071878.1.gb-7-CDS]
gx/coords/test_coords_export.py::test_co_5_exported_location_string_still_carries_the_angle_bracket[NC_001699.1.gb-21-exon]
gx/coords/test_coords_export.py::test_co_6_order_is_not_rewritten_as_join
gx/coords/test_coords_export.py::test_co_7_remote_reference_is_not_silently_localised
```

## Acceptance

- The twelve sealed gx tests listed above pass, unedited.
- `make check` stays green.
- Nothing is removed. In particular the currently-passing requirements over the same code path keep passing: CO-2 (export is the identity on every interval, part order, strand and extracted subsequence for both real records), CO-3 (`join()` parts survive load and export, including the origin-spanning `join(5118..5130,1..12)` and `join(1975..1993,1..19)`, with `circular` topology preserved), CO-4 (`complement` and `complement(join)` keep strand and part order), and the `join(101..110,141..150)` control in the CO-6 probe — a fix for `order()` must not be made by breaking `join()`.
