# Zero branch length clamping

[`fix_branch_length()`](../../packages/treetime/src/hacks/fix_branch_length.rs#L6) (`#fix_branch_length`) clamps branch lengths to `MIN_BRANCH_LENGTH_FRACTION * (1/L)` to avoid NaN in the marginal backward pass, which computes `log(expQt)` -- divergent at `t=0`. Mathematically, `expQt(0) = I` is correct (identity matrix), but the log-space computation path requires `t > 0`.

This is a documented limitation matching v0 behavior. The consequence: [test assertions for zero-branch-length cases](../../packages/treetime/src/ancestral/__tests__/test_marginal_dense.rs) were weakened from exact mathematical properties (`Ti == 0` with epsilon 1e-15, `off_diagonal < 0.1`) to `is_finite()` checks.

## Related issues

- Source: [N-core-branch-length-clamping.md](../issues/N-core-branch-length-clamping.md) -- delete after full resolution
