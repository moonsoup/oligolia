# gx trace — restriction/digest, Cycle A round 1

`gx/tests.json` carries the contract's `gx-re-<n>` ids; pytest reports node ids.
This table is the mapping between them, in declaration order per file. Test
function names carry the *requirement* id (`test_re_<n>_…`), so a failing node id
names its requirement directly; the `gx-re-<n>` column is the register-side id.

| gx id | node | covers |
|---|---|---|
| gx-re-1 | `gx/test_restriction_sites.py::test_re_1_site_starts_on_linear_puc19` | RE-1 |
| gx-re-2 | `gx/test_restriction_sites.py::test_re_1_count_and_pattern_fields` | RE-1 |
| gx-re-3 | `gx/test_restriction_sites.py::test_re_2_panel_recognition_sites_are_palindromic` | RE-2 |
| gx-re-4 | `gx/test_restriction_sites.py::test_re_2_reverse_strand_gives_mirror_positions` | RE-2 |
| gx-re-5 | `gx/test_restriction_sites.py::test_re_3_site_on_bottom_strand_is_found` | RE-3 |
| gx-re-6 | `gx/test_restriction_sites.py::test_re_3_site_count_matches_cut_count` | RE-3, RE-12 |
| gx-re-7 | `gx/test_restriction_sites.py::test_re_10_ambiguous_avai_all_four_resolutions` | RE-10 |
| gx-re-8 | `gx/test_restriction_sites.py::test_re_10_ambiguous_avaii_both_resolutions` | RE-10 |
| gx-re-9 | `gx/test_restriction_sites.py::test_re_10_ambiguous_enzymes_on_puc19` | RE-10 |
| gx-re-10 | `gx/test_restriction_cut_convention.py::test_re_4_cut_position_convention` (5 params) | RE-4 |
| gx-re-11 | `gx/test_restriction_cut_convention.py::test_re_4_whole_panel_cut_positions` | RE-4 |
| gx-re-12 | `gx/test_restriction_cut_convention.py::test_re_9_overlapping_sites_are_each_reported` | RE-9 |
| gx-re-13 | `gx/test_restriction_cut_convention.py::test_re_9_overlapping_sites_of_different_enzymes` | RE-9 |
| gx-re-14 | `gx/test_restriction_cut_convention.py::test_re_14_digest_is_independent_of_enzyme_order` | RE-14 |
| gx-re-15 | `gx/test_restriction_digest_fragments.py::test_re_6_linear_fragments_match_biopython` (3 params) | RE-6 |
| gx-re-16 | `gx/test_restriction_digest_fragments.py::test_re_6_linear_whole_panel` | RE-6 |
| gx-re-17 | `gx/test_restriction_digest_fragments.py::test_re_7_circular_fragments_match_biopython` (3 params) | RE-7 |
| gx-re-18 | `gx/test_restriction_digest_fragments.py::test_re_7_circular_whole_panel` | RE-7 |
| gx-re-19 | `gx/test_restriction_digest_fragments.py::test_re_8_enzyme_with_no_site` (20 params) | RE-8 |
| gx-re-20 | `gx/test_restriction_circular_origin.py::test_re_5_origin_spanning_site_found_exactly_once` | RE-5 |
| gx-re-21 | `gx/test_restriction_circular_origin.py::test_re_5_origin_spanning_site_absent_when_linear` | RE-5 |
| gx-re-22 | `gx/test_restriction_circular_origin.py::test_re_5_origin_spanning_cut_and_single_fragment` | RE-5, RE-7 |
| gx-re-23 | `gx/test_restriction_circular_origin.py::test_re_5_origin_spanning_site_found_once_for_every_rotation` | RE-5 |
| gx-re-24 | `gx/test_restriction_input_normalisation.py::test_re_11_layout_characters_do_not_move_the_coordinates` (3 params) | RE-11 |
| gx-re-25 | `gx/test_restriction_input_normalisation.py::test_re_11_fragment_sequences_are_nucleotides_only` (3 params) | RE-11 |
| gx-re-26 | `gx/test_restriction_input_normalisation.py::test_re_12_site_and_cut_endpoints_share_one_coordinate_system` (3 params) | RE-12 |
| gx-re-27 | `gx/test_restriction_input_normalisation.py::test_re_12_coherence_holds_across_the_panel` | RE-12 |
| gx-re-28 | `gx/test_restriction_input_normalisation.py::test_re_13_non_iupac_template_is_rejected_by_both_endpoints` | RE-13 |
| gx-re-29 | `gx/test_restriction_input_normalisation.py::test_re_13_unknown_enzyme_is_the_reference_rejection` | RE-13 |

## Running

```bash
cd /workspace/oligolia-gx && PYTHONPATH=. backend/.venv/bin/python -m pytest gx/ -q -p no:cacheprovider
```

Offline: the only sequence input is `gx/corpus/inputs/puc19_L09137.2.txt`; the
oracles are `Bio.Restriction` (pinned in `gx/corpus/`) and values re-derived by
index. `gx/` is outside `backend/tests/` and `gui/`, so `make check` — which
lints `backend/ gui/ oligolia.py .claude/qa/` and runs those two suites — is
unaffected.
