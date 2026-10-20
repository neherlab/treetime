# Test quality deficiencies

## Summary

Remaining test quality issues: unsourced and symmetric-model tests, and weak assertions.

## Instances

### propagate_raw_per_site tests use a symmetric model and unsourced expected values

[packages/treetime/src/partition/marginal/sparse/__tests__/test_message.rs#L17](../../packages/treetime/src/partition/marginal/sparse/__tests__/test_message.rs#L17), [#L73](../../packages/treetime/src/partition/marginal/sparse/__tests__/test_message.rs#L73)

Both test functions use JC69 with per-site rates, whose transition matrix is symmetric ($P(t) = P(t)^T$). The forward (`transpose = false`) and backward (`transpose = true`) tests therefore produce identical expected values and cannot detect a defect in the `transpose` handling of `fn propagate_raw_per_site()` [packages/treetime/src/partition/marginal/sparse/message.rs#L137](../../packages/treetime/src/partition/marginal/sparse/message.rs#L137). The expected values are hardcoded literals with no stated source ([#L39-L56](../../packages/treetime/src/partition/marginal/sparse/__tests__/test_message.rs#L39-L56), [#L95-L112](../../packages/treetime/src/partition/marginal/sparse/__tests__/test_message.rs#L95-L112)), so whether they were captured from the system under test is unknown.

Fix:

- Add an HKY85 (or other non-symmetric model) case where forward and backward produce different results, so the transpose code path is exercised. HKY85 with unequal base frequencies (e.g. $\pi = [0.4, 0.1, 0.1, 0.4]$) produces an asymmetric $P(t)$
- Derive expected values from the HKY85 closed form or the v0 Python GTR, never from the v1 GTR matrix exponential, because that reintroduces a circular oracle
- State the source of the expected values in the test

### Skyline optimizer tests check structure only

[packages/treetime/src/coalescent/__tests__/test_skyline.rs](../../packages/treetime/src/coalescent/__tests__/test_skyline.rs)

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
