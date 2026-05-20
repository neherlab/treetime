# Test quality deficiencies

## Summary

Remaining test quality issues: circular tests, weak assertions, and missing coverage for specific entities.

## Instances

### propagate_raw_per_site tests are circular

`packages/treetime/src/partition/marginal_helpers.rs:209,245:`

Both test functions compute expected values using `gtr.expQt_with_rate()`, the same function called internally by the SUT. Tautological verification.

### Skyline tests are runs-to-completion only

`packages/treetime/src/coalescent/__tests__/test_skyline.rs`

Assert finite/positive only, no numerical verification against known values.

### Polytomy topology tests use point distributions

`packages/treetime/src/timetree/optimization/__tests__/test_polytomy.rs:46:`

Branch-length distributions are point masses making the cost function trivial: any merge time not exactly at the point yields `ln(1e-10)` from the `unwrap_or` fallback.

### test_gm_runner_marginal_sparse compares against dense expected (NOT A DEFECT)

Cross-mode validation: v0 has no sparse mode, so dense oracle is the only available comparison. Documented in test support module.

### Missing test coverage for specific entities

- No tests for `fn Sub::from_str`, `fn parse_pos`, validators at `seq/mutation.rs`
- No tests for `enum AlphabetName::AaNoStop` at `alphabet.rs`
