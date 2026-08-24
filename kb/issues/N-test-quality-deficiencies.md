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

### Missing test coverage for specific entities

- No tests for `fn Sub::from_str`, `fn parse_pos`, validators at `seq/mutation.rs`
- No tests for `enum AlphabetName::AaNoStop` at `alphabet.rs`

## Related tickets

- [kb/tickets/test-add-hky85-case-to-propagate-raw-per-site.md](../tickets/test-add-hky85-case-to-propagate-raw-per-site.md)
