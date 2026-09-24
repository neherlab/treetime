# Test quality deficiencies

## Summary

Remaining test quality issues: circular tests, weak assertions, and missing coverage for specific entities.

## Instances

### propagate_raw_per_site tests are circular

`packages/treetime/src/partition/marginal_helpers/__tests__/test_marginal_helpers.rs:10,66:`

Both test functions compute expected values using `gtr.expQt_with_rate()`, the same function called internally by the SUT. Tautological verification.

### Skyline tests are runs-to-completion only

`packages/treetime/src/coalescent/__tests__/test_skyline.rs`

Assert finite/positive only, no numerical verification against known values.

### test_gm_runner_marginal_sparse compares against dense expected (NOT A DEFECT)

Cross-mode validation: v0 has no sparse mode, so dense oracle is the only available comparison. Documented in test support module.

### Jukes-Cantor distance tests use loose bounds and unsourced values

`packages/treetime/src/gtr/__tests__/test_jc_distance.rs`

- `test_jukes_cantor_distance_saturation_cap_order_of_magnitude` and `test_jukes_cantor_distance_issue_documented_error` accept ranges such as `10.0..11.0` and `0.06..0.08` where the closed form $d = -\frac{k-1}{k} \ln\left(1 - \frac{k}{k-1} p\right)$ gives exact values
- `test_jukes_cantor_distance_small_p_approaches_p` accepts a relative error below `1e-3`; the first-order expansion $d \approx p + \frac{k}{2(k-1)} p^2$ supports a much tighter bound
- `test_jukes_cantor_distance_known_values` compares at `1e-15` against values with no stated source
- The monotonicity and `d >= p` tests loop over sampled inputs with assertions inside the loop; a property test states the invariant directly

### Poisson indel tests check signs instead of values

`packages/treetime/src/optimize/__tests__/test_indel.rs`

- `test_optimize_indel_poisson_derivative_positive_near_zero` asserts `derivative > 1e5`; the exact derivative $k/t - \mu$ is $999995$ for its inputs
- `test_optimize_indel_poisson_second_derivative_negative` asserts only the sign of $-k/t^2$
- `test_optimize_indel_statrs_ln_factorial` tests the third-party `statrs::function::factorial::ln_factorial` instead of project code

### Error tests match message substrings

These tests assert `err.to_string().contains(...)` instead of the exact message with `assert_error!`, so a changed message or a different error with the same fragment still passes:

- `packages/app-cli/src/commands/ancestral/__tests__/test_aa_node_data.rs`: `test_validate_aa_args_requires_cds_placeholder`, `test_validate_aa_args_empty_cdses_no_annotation_errors`, `test_validate_aa_root_sequence_cdses_requires_every_cds`
- `packages/treetime-io/src/__tests__/test_gff.rs`: `test_parse_gff3_cds_features_rejects_non_multiple_of_three`, `test_parse_gff3_cds_features_rejects_mixed_seqids`
- `packages/treetime-graph/src/__tests__/test_topology_order.rs`: `topology_order_rejects_cycles`, `topology_order_target_order_rejects_empty`, `topology_order_target_order_rejects_duplicate_ranking_labels`, `topology_order_target_order_rejects_duplicate_final_leaf_labels`

### Missing test coverage for specific entities

- No tests for `fn Sub::from_str`, `fn parse_pos`, validators at `seq/mutation.rs`
- No tests for `enum AlphabetName::AaNoStop` at `alphabet.rs`

## Related tickets

- [kb/tickets/test-add-hky85-case-to-propagate-raw-per-site.md](../tickets/test-add-hky85-case-to-propagate-raw-per-site.md)
