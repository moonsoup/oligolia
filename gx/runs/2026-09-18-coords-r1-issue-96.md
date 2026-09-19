The HTTP edit endpoint builds its stored record from scratch and passes only five fields. `annotations` and `is_circular` are not among them, so both fall back to their model defaults: an empty list and `False`. Upload an annotated circular plasmid, insert one base through the API, download the GenBank — every annotation is gone and the molecule is linear, with no error and no warning.

Found by the Rockin-Robin `oligolia-gx-tests` contract, coordinates round 1, against commit `4bb77b3`.

## Requirement

**CO-14** — *The sequence-edit endpoint must not silently discard the record's annotations and topology. `POST /sequences/{id}/edit` stores the edited molecule under a new id; that stored record must carry the record's features — shifted for the edit, the way the GUI's own edit path does via `shift_annotations` — and must keep `is_circular`, because the stored record is what the file endpoints then export.*

## Root cause

`backend/routers/sequences.py:135-145`:

```python
# Persist the edited sequence under a new ID
new_id = f"{seq_id}_{op}"
_store[new_id] = Sequence(
    id=new_id,
    name=f"{s.name} ({op})",
    description=s.description,
    seq=result.result_seq,
    molecule_type=MoleculeType.PROTEIN if op == "translate" else
                  MoleculeType.RNA if op == "transcribe" else
                  MoleculeType.DNA if op == "back_transcribe" else s.molecule_type,
)
```

Five fields. `annotations` and `is_circular` are never passed, so `annotations=[]` and `is_circular=False`. Neither `shift_annotations` nor `flip_annotations` is called anywhere in `backend/routers/` — the only call sites are `gui/panels/sequence_panel.py:959` and `:964`.

## Differential evidence

Input: `gx/corpus/inputs/NC_001699.1.gb` (JC polyomavirus, 5130 bp, **circular**, 28 features), sha256 `6935cc5c717c934b1504a13d8a4c2f457c80adfbe7e9d0c16f1971a2d67a7dd0`, loaded and stored under id `NC_001699.1`.

| operation | expected stored record | actual stored record |
|---|---|---|
| `insert` at 3000, 30 × `A` → `NC_001699.1_insert` | 5160 bp, 28 features shifted, `is_circular=True` | 5160 bp, **0 features**, `is_circular=False` |
| `delete` 2100–2200 → `NC_001699.1_delete` | 5030 bp, features shifted | 5030 bp, **0 features**, `is_circular=False` |
| `replace` 1600–1650 with 10 × `A` → `NC_001699.1_replace` | 5090 bp, features shifted | 5090 bp, **0 features**, `is_circular=False` |
| `reverse_complement` → `NC_001699.1_reverse_complement` | 5130 bp, 28 features flipped | 5130 bp, **0 features**, `is_circular=False` |

Literally: `Sequence(id='NC_001699.1_insert', …, annotations=[], …)`. 0 of 28, four operations, four identical results — including `reverse_complement`, the one case where the GUI does have `flip_annotations`.

The **bases** the endpoint returns are correct for all four operations (the control assertion in the same file passes), so this is about the stored record, not the splice.

End of the chain: `write_genbank` of the stored record produces an **empty FEATURES block**. The stored record is what `GET /sequences/{id}` returns and what `POST /files/download/genbank` writes, so edit-then-export over the API is a documented way to strip a GenBank file of its annotations.

Topology is not cosmetic here. The repo's own restriction/digest work keys off `is_circular` — an origin-spanning cut only exists on a circular molecule — and the Sequence tab's topology toggle is documented as *"preserved across edits and GenBank export"* (`gui/panels/sequence_panel.py:466`). One insert through the API silently linearises the plasmid.

## Relationship to closed issues

This is **#59 left unfixed on the HTTP side**, and it is the wider half. #59 ("Sequence edits never shift annotations; … annotations are never shifted. Insert, delete, replace and reverse-complement all leave feature coordinates untouched") was fixed by giving the GUI path `shift_annotations`/`flip_annotations`; the router path was never given anything. Where #59's GUI symptom was features pointing at the wrong bases, the router symptom is that there are no features at all — it removes every annotation rather than misplacing one.

## Pinned source anchors

| anchor | sha256 |
|---|---|
| `biopython-1.85/Bio/SeqFeature.py#    def extract(self, parent_sequence, references=None):` | `9bdb8d0386de1848092e8c7901e80b0955074ff0cafc7e89c4e41ae4d07b7d56` |
| `inputs/NC_001699.1.gb#     CDS             complement(join(2603..4426,4771..5013))` | `6935cc5c717c934b1504a13d8a4c2f457c80adfbe7e9d0c16f1971a2d67a7dd0` |

## Sealed failing tests

Seal digest: `2ef985275ccf556d2f29023f030d029fabb0a5e5d723ad3c33591a33765987f9` (`gx/coords/seal.json`). Same 42 failed / 158 passed on Linux and macOS. **The sealed tests must not be edited.**

```
gx/coords/test_coords_edit_endpoint.py::test_co_14_the_stored_record_keeps_its_annotations[insert]
gx/coords/test_coords_edit_endpoint.py::test_co_14_the_stored_record_keeps_its_annotations[delete]
gx/coords/test_coords_edit_endpoint.py::test_co_14_the_stored_record_keeps_its_annotations[replace]
gx/coords/test_coords_edit_endpoint.py::test_co_14_the_stored_record_keeps_its_annotations[reverse_complement]
gx/coords/test_coords_edit_endpoint.py::test_co_14_the_stored_record_keeps_its_topology[insert]
gx/coords/test_coords_edit_endpoint.py::test_co_14_the_stored_record_keeps_its_topology[delete]
gx/coords/test_coords_edit_endpoint.py::test_co_14_the_stored_record_keeps_its_topology[replace]
gx/coords/test_coords_edit_endpoint.py::test_co_14_the_stored_record_keeps_its_topology[reverse_complement]
gx/coords/test_coords_edit_endpoint.py::test_co_14_surviving_annotations_point_at_the_right_bases
gx/coords/test_coords_edit_endpoint.py::test_co_14_a_downloaded_genbank_of_the_edited_record_still_has_features
```

Note on `test_co_14_surviving_annotations_point_at_the_right_bases`: today it fails because nothing survives, so it only restates the loss. Once the endpoint carries annotations, it becomes the check that they were *shifted* — it asks whether a surviving feature still extracts bases the original record contained.

## Acceptance

- The ten sealed gx tests listed above pass, unedited.
- `make check` stays green, including the existing `/sequences/{id}/edit` coverage in the backend suite.
- Nothing is removed: the endpoint's returned bases stay correct for all four operations (the control that passes today must keep passing), the `translate`/`transcribe`/`back_transcribe` molecule-type handling at `backend/routers/sequences.py:142-144` is unchanged, and the new-id scheme (`{seq_id}_{op}`) is unchanged.
