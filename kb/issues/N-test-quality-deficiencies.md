# Test quality deficiencies

## Summary

Remaining test quality issues: unsourced expected values and weak assertions.

## Instances

### Skyline optimizer tests check structure only

[`packages/treetime/src/coalescent/__tests__/test_skyline.rs`](../../packages/treetime/src/coalescent/__tests__/test_skyline.rs)

Two star-tree tests check exact values: `test_skyline_constant_tc_of_a_star_tree_is_the_kingman_estimate` compares the constant Tc and its log-likelihood with the closed-form Kingman estimate, and `test_skyline_merger_rate_and_span_of_a_star_tree_follow_the_lineage_count` checks the merger rate against the lineage count. The tests of `optimize_skyline()` assert only segment counts, finite log-likelihoods, positive Tc values, the comparison with constant Tc, and scale invariance. No test compares the multi-segment Tc trajectory with known values.

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

### Maximum grid point tests compare error text through `unwrap_err()`

- [`packages/treetime-grid/src/__tests__/max_grid_points.rs#L39-L40`](../../packages/treetime-grid/src/__tests__/max_grid_points.rs#L39-L40): `test_max_grid_points_from_str_rejects` compares `MaxGridPoints::from_str(text).unwrap_err().to_string()` with the expected message
- [packages/app-server/src/bin/treetime-server.rs#L190-L192](../../packages/app-server/src/bin/treetime-server.rs#L190-L192): `test_treetime_server_max_grid_points_out_of_range_fails` compares `ServerArgs::try_parse_from(..).unwrap_err().to_string()` with the expected message

Both tests panic without a message when the call succeeds. `assert_error!` of `packages/treetime-utils/src/testing/assert.rs` checks the `Err` and the exact message in one assertion and reports an unexpected `Ok`.
