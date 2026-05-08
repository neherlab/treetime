# Test quality deficiencies

## Summary

Systematic test quality issues: missing pretty_assertions, missing rstest trace attributes, mixed assertion macro styles, loop-over-cases anti-pattern, helper placement violations, circular tests, weak assertions, loose tolerances, and other anti-patterns.

## Instances

### Missing pretty_assertions (RESOLVED)

### Missing #[trace] on rstest (RESOLVED)

### Mixed approx wrappers (RESOLVED)

### Loop-over-cases in test assertions (RESOLVED)

Converted element-by-element zip assertions to whole-array macros (33 instances across 6 files) and loop-over-test-cases to rstest parameterized tests (10 instances across 4 files). Remaining loops in proptest files are structural invariant checks with per-element diagnostics, not anti-patterns.

### Sparse root-invariance proptest tolerance (RESOLVED)

Tightened to 1e-6 and `#[ignore]`d pending fix of `M-ancestral-sparse-root-invariance.md`.

### propagate_raw_per_site tests are circular

`packages/treetime/src/representation/partition/marginal_helpers.rs:209,245:`

Both test functions compute expected values using `gtr.expQt_with_rate()`, the same function called internally by the SUT. Tautological verification.

### Element-by-element float loops in alphabet tests (RESOLVED)

Replaced 22 element-by-element loops with `pretty_assert_ulps_eq!` on whole arrays.

### GM runner tests use grossly loose tolerance (RESOLVED)

Both tests tightened to `epsilon = 1e-6` and `#[ignore]`d pending fix of `kb/issues/M-timetree-branch-grid-uniform-resolution.md`.

### Skyline tests are runs-to-completion only

`packages/treetime/src/commands/timetree/coalescent/__tests__/test_skyline.rs`

Assert finite/positive only, no numerical verification against known values.

### Polytomy topology tests use point distributions

`packages/treetime/src/commands/timetree/optimization/__tests__/test_polytomy.rs:46:`

Branch-length distributions are point masses making the cost function trivial: any merge time not exactly at the point yields `ln(1e-10)` from the `unwrap_or` fallback.

### test_gm_runner_marginal_sparse compares against dense expected (NOT A DEFECT)

Cross-mode validation: v0 has no sparse mode, so dense oracle is the only available comparison. Documented in test support module.

### Loose tolerance 1e-6 in div tests (RESOLVED)

Tightened to measured values: 5 tests at `1e-8`, deep_tree at `1e-7`.

### Missing test coverage for specific entities

- No tests for `fn Sub::from_str`, `fn parse_pos`, validators at `seq/mutation.rs`
- No tests for `enum AlphabetName::AaNoStop` at `alphabet.rs`

## Tolerance violations

### Grossly loose (RESOLVED)

- ~~`commands/timetree/coalescent/__tests__/test_integration.rs:157:` `epsilon = 10.0`~~ Tightened to 1e-6, ignored (midpoint-rule discretization)
- ~~`commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs:80:` `epsilon = 1e0`~~ Tightened to 1e-6, ignored

### Non-standard format (RESOLVED)

- ~~`commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs:46:` `epsilon = 3e-1`~~ Tightened to 1e-6, ignored (grid resolution)
- ~~`commands/mugration/__tests__/test_gm_mugration.rs:92:` `epsilon = 2e-2`~~ Tightened to 1e-6, ignored (D1/D2 divergence)

### Variable tolerances (RESOLVED)

All replaced with strict numeric literals:

- ~~`GRID_SPACING_BL`, `expected_epsilon`~~: now `1e-10`, `1e-2` (measured)
- ~~`epsilon = dx`~~: now `1e-3` (measured)
- ~~`PASS_THRESHOLD`~~: now inline `1e-5`
- ~~`epsilon = 1e-6` variable~~: now `max_ulps = 4` (measured 1 ULP)
- ~~`1e-2` with `TODO(investigate)`~~: kept 1e-2 (grid accuracy limit), removed TODO
