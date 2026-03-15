# Standalone sparse optimizer skips derivative check for zero branches

`run_optimize_sparse()` ([packages/treetime/src/commands/optimize/optimize_sparse.rs#L142-L222](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L142-L222)) uses a stricter zero-branch likelihood threshold (0.0001) than the unified path (0.01 in `is_zero_branch_optimal()`) and does not check the derivative sign before setting branch length to zero. A TODO comment at [optimize_sparse.rs#L172](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L172) acknowledges this: "could check that derivative is negative."

The unified path in `is_zero_branch_optimal()` ([packages/treetime/src/commands/optimize/optimize_unified.rs#L153-L173](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L153-L173)) checks both conditions: likelihood at zero > 0.01 AND derivative < 0. The standalone sparse function short-circuits to zero when likelihood > 0.0001 alone, which can set branches to zero even when the likelihood surface slopes upward at zero (meaning a positive branch length would be better).

`run_optimize_sparse()` and its dense counterpart `run_optimize_dense()` are not called by the main `optimize` command. They serve as test-only helpers (called from `__tests__/test_dense_sparse_equivalence/`). The inconsistency is low-risk but makes equivalence tests unreliable as a correctness check for the unified path.

## Proposed solution

Either align `run_optimize_sparse()` with the unified path (add derivative check, raise threshold to 0.01) or mark both standalone functions as `#[cfg(test)]` since they only serve tests.
