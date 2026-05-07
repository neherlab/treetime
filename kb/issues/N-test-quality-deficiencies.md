# Test quality deficiencies

## Summary

Systematic test quality issues: missing pretty_assertions, missing rstest trace attributes, mixed assertion macro styles, loop-over-cases anti-pattern, helper placement violations, circular tests, weak assertions, loose tolerances, and other anti-patterns.

## Instances

### Missing pretty_assertions (8 files, 46+ raw assert_eq!)

1. `commands/ancestral/__tests__/test_fitch.rs` (14)
2. `commands/clock/__tests__/test_date_constraints.rs` (14)
3. `commands/clock/__tests__/test_clock_filter.rs` (2)
4. `commands/clock/__tests__/test_clock_dengue100.rs` (1)
5. `commands/ancestral/__tests__/test_marginal_consistency.rs` (3)
6. `seq/__tests__/test_find_char_ranges.rs` (4)
7. `seq/__tests__/test_div.rs` (7)
8. `commands/mugration/__tests__/test_comment_output.rs` (1)

### Missing #[trace] on rstest (5 files, 17 blocks)

1. `gtr/__tests__/test_gm_gtr.rs` (7)
2. `gtr/__tests__/test_gm_gtr_site_specific.rs` (3)
3. `gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_parameterized.rs` (2)
4. `representation/__tests__/test_partition_marginal_sparse.rs` (1)
5. `seq/__tests__/test_find_char_ranges.rs` (4)

### Mixed approx wrappers (4 files)

Files using both `approx::assert_ulps_eq!` and `pretty_assert_ulps_eq!`:

1. `commands/ancestral/__tests__/test_marginal_dense.rs`
2. `commands/ancestral/__tests__/test_marginal_sparse.rs`
3. `commands/ancestral/__tests__/test_marginal_consistency.rs`
4. `commands/clock/__tests__/test_clock_regression.rs`

### Loop-over-cases in test assertions (21 files)

Element-by-element loop assertions. Highest density: `test_prop_gtr_site_specific.rs` (39 loops), `generators.rs` (13 loops). Many within proptest bodies where `prop_assert!` is correct but `prop_assert_array_*_eq!` from project utilities could replace the element iteration.

### Helpers placement violations (9 files)

Helper functions before tests and/or not wrapped in `mod helpers`:

1. `seq/__tests__/test_mutation.rs`
2. `seq/__tests__/test_find_char_ranges.rs`
3. `alphabet/__tests__/test_alphabet_config.rs`
4. `commands/timetree/coalescent/__tests__/test_skyline.rs`
5. `commands/timetree/optimization/__tests__/test_relaxed_clock.rs`
6. `commands/timetree/optimization/__tests__/test_polytomy.rs`
7. `commands/timetree/convergence/__tests__/test_metrics.rs`
8. `commands/timetree/inference/__tests__/test_branch_length_likelihood.rs`
9. `gtr/__tests__/test_gm_gtr_site_specific.rs:170-207:`

### Sparse root-invariance proptest tolerance 1e-1

`packages/treetime/src/commands/ancestral/__tests__/test_marginal_root_invariance_prop.rs:57:`

Hides >2-orders-of-magnitude pulley-principle violation. Dense uses 1e-6. The 5-order gap indicates real algorithmic divergence. Related: `M-ancestral-sparse-root-invariance.md`.

### propagate_raw_per_site tests are circular

`packages/treetime/src/representation/partition/marginal_helpers.rs:209,245:`

Both test functions compute expected values using `gtr.expQt_with_rate()`, the same function called internally by the SUT. Tautological verification.

### Element-by-element float loops in alphabet tests

`packages/treetime/src/alphabet/__tests__/` - 19+ instances of `for (actual, expected) in ... { assert_ulps_eq!(...) }`. Should use `pretty_assert_ulps_eq!` or `pretty_assert_abs_diff_eq!` on whole arrays.

### GM runner tests use grossly loose tolerance (RESOLVED)

Both tests tightened to `epsilon = 1e-6` and `#[ignore]`d pending fix of
`kb/issues/M-timetree-branch-grid-uniform-resolution.md`.

### Skyline tests are runs-to-completion only

`packages/treetime/src/commands/timetree/coalescent/__tests__/test_skyline.rs`

Assert finite/positive only, no numerical verification against known values.

### Polytomy topology tests use point distributions

`packages/treetime/src/commands/timetree/optimization/__tests__/test_polytomy.rs:46:`

Branch-length distributions are point masses making the cost function trivial: any merge time not exactly at the point yields `ln(1e-10)` from the `unwrap_or` fallback.

### test_gm_runner_marginal_sparse compares against dense expected (NOT A DEFECT)

Cross-mode validation: v0 has no sparse mode, so dense oracle is the only available comparison. Documented in test support module.

### .read() instead of .read_arc() (8 instances)

`packages/treetime/src/commands/ancestral/__tests__/test_python_parity.rs:192,240,283,327,390,391,463,473:`

Test code calls `.read()` on `parking_lot::RwLock` values wrapped in `Arc`. Per project convention, `.read_arc()` should be used to return an `ArcRwLockReadGuard` that keeps the `Arc` alive.

### Loose tolerance 1e-6 in div tests

`packages/treetime/src/seq/__tests__/test_div.rs`

6 assertions use `epsilon = 1e-6`, 1 uses `epsilon = 1e-9`. Whether 1e-6 is the tightest passing tolerance has not been measured.

### Missing test coverage for specific entities

- No tests for `fn Sub::from_str`, `fn parse_pos`, validators at `seq/mutation.rs`
- No tests for `enum AlphabetName::AaNoStop` at `alphabet.rs`

## Tolerance violations

### Grossly loose (hard failure per project rules)

- `commands/timetree/coalescent/__tests__/test_integration.rs:157:` `epsilon = 10.0`
- ~~`commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs:80:` `epsilon = 1e0`~~ RESOLVED: tightened to 1e-6, ignored

### Non-standard format

- `commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs:46:` `epsilon = 3e-1`
- `commands/mugration/__tests__/test_gm_mugration.rs:92:` `epsilon = 2e-2`

### Variable tolerances

- `commands/timetree/inference/__tests__/test_branch_length_likelihood.rs:80,129,162:` `GRID_SPACING_BL`, `expected_epsilon` (computed)
- `commands/timetree/output/__tests__/test_confidence_extract.rs:232:` `epsilon = dx` (computed)
- `commands/timetree/coalescent/__tests__/test_gm_coalescent.rs:102:` `worst_err < PASS_THRESHOLD` (named constant)
- `commands/ancestral/__tests__/test_marginal_dense.rs:364:` `let epsilon = 1e-6;` (variable)
- `gtr/__tests__/test_prop_gtr_site_specific.rs:283:` `1e-2` with `TODO(investigate)` tag. Justification: linear interpolation on 61-point grid has inherent accuracy limits; observed max error ~1.8e-3.
