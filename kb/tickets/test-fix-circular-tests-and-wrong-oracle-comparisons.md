# Fix circular tests and wrong oracle comparisons

## Circular propagate_raw_per_site tests

`packages/treetime/src/representation/partition/marginal_helpers.rs:209,245:`

Both test functions compute expected values using `gtr.expQt_with_rate()`, the same function called internally by the SUT. Tautological verification.

## test_gm_runner_marginal_sparse compares against dense expected

`packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs:45:`

Cross-mode validation (no independent sparse oracle exists).

## Related issues

- Source: [N-test-quality-deficiencies.md](../issues/N-test-quality-deficiencies.md) -- delete after full resolution
