# Match v0 coarse tolerance in optimize pre-step

v0's `optimize_tree_marginal` uses a progressive tolerance schedule: `tol = 1e-8 + 0.01^(i+1)` where `i` is the iteration index (0-based). At iteration 0 (the pre-step), this evaluates to `tol = 0.02`, intentionally coarse to prioritize speed over precision in early iterations.

v1's pre-step uses `BrentSqrt` with its default tolerance (Brent's method default from `argmin`, approximately 1e-8). This is tighter than v0's first-iteration tolerance.

## Impact

Low. The tighter tolerance may cause the Brent optimizer to take a few extra function evaluations per edge during the pre-step, slightly increasing runtime. Branch length values may differ from v0 at the level of the tolerance difference (0.02 vs 1e-8), but this is within the noise of subsequent time inference refinement. Golden master comparison may show small differences in pre-step branch lengths.

## Proposed solution

Add a `tolerance` parameter to `run_optimize_mixed()` or pass it through `BranchOptMethod`. The pre-step would use `tol=0.02` (matching v0 iteration 0), while the optimize command's loop would use the progressive schedule.

Requires plumbing the tolerance through the Brent method implementation in `method_brent.rs`, which currently uses `argmin::BrentOpt` with hardcoded tolerance.

## Related issues

- Source: [kb/issues/N-optimize-pre-step-coarse-tolerance.md](../issues/N-optimize-pre-step-coarse-tolerance.md) -- delete after full resolution
