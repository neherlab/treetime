# Per-edge optimization: dispatch zero boundary check

## Problem

All 6 inner-loop functions (Brent, BrentSqrt, BrentLog, Newton, NewtonSqrt, NewtonLog) are wired in `run_optimize_mixed`. One cross-cutting dispatch-level concern remains:

**Zero boundary gap**: methods that cannot reach $t = 0$ internally (NewtonLog domain excludes zero; Brent variants use a strictly positive lower bracket) miss the boundary optimum on non-unimodal models where the `is_zero_branch_optimal()` shortcut is skipped.

## Background

The `is_zero_branch_optimal()` function is a fast pre-dispatch check that returns true when $\ell'(0) < 0$. It is restricted to models with unimodal branch-length likelihood: JC69, F81, and binary models (Dinh & Matsen 2017, Corollary 3.1). For these models, a negative derivative at $t=0$ proves zero is the global maximum because no other stationary point can exist on $(0, \infty)$.

For non-unimodal models (K80, HKY85, TN93, general GTR), multiple local maxima can exist, so the derivative-sign criterion is insufficient. The shortcut returns false and dispatch falls through to the chosen optimizer. The optimizer then runs on the open domain $t > 0$.

Methods that cannot evaluate exactly at $t = 0$:

- **NewtonLog**: $u = \ln(t)$ is undefined at zero. The implementation bumps zero starts to `one_mutation` and uses $u_{\min} = \ln(\max(\text{min\_bl}, 10^{-12}))$ as a finite floor.
- **Brent variants**: bracket lower bound is `max(min_branch_length, 1e-12)`, strictly positive.

For non-unimodal models, edges whose true optimum is exactly $t = 0$ are biased to a tiny positive value. Downstream effect: `find_zero_optimal_internal_edges()` (in `run.rs`) never collects them for collapse, so the optimize loop's topology cleanup misses these edges. This changes default optimize behavior for non-unimodal models.

The grid search fallback (`grid_search_inner`) does compare zero against the grid best via `is_zero_better_than_grid_best()`, but the Newton and Brent paths have no equivalent post-optimization check.

## Approach

Add a post-optimization zero-comparison in `run_optimize_mixed`, mirroring the grid search boundary check. After any method returns $t^*$:

```rust
let new_branch_length = match method { ... };

// Post-optimization boundary check: compare against t=0 for methods
// that cannot reach zero internally (NewtonLog, Brent variants).
let new_branch_length = if new_branch_length > 0.0
  && is_zero_better_than_grid_best(contributions, indel_count, indel_rate, new_branch_length)
{
  0.0
} else {
  new_branch_length
};
```

This applies to all methods uniformly. For Newton and NewtonSqrt (which can reach zero through their clamps), it is a no-op when the optimizer already found zero. For NewtonLog and the Brent variants, it provides the missing boundary comparison.

The check is gated on `indel_count == 0` and `all_sites_valid_at_zero` (handled inside `is_zero_better_than_grid_best`). When indels are present, the Poisson log-likelihood at $t = 0$ is $-\infty$ so zero is never optimal.

## Verification

`./dev/docker/run ./dev/dev t` -- all existing tests pass.

Boundary check regression test: configure a non-unimodal GTR model where $t = 0$ is the true optimum on a chosen edge. Verify all 6 methods return zero (not a tiny positive value).

## Dependencies

- Depends on: nothing further (all inner loops are wired).
- Depended on by: [M-optimize-method-tests](M-optimize-method-tests.md), [L-optimize-method-docs](L-optimize-method-docs.md)

## Cross-references

- Dispatch: `packages/treetime/src/commands/optimize/optimize_unified.rs` (`run_optimize_mixed`)
- Zero-branch shortcut: `packages/treetime/src/commands/optimize/optimize_unified.rs` (`is_zero_branch_optimal`)
- Grid zero comparison: `packages/treetime/src/commands/optimize/optimize_unified.rs` (`is_zero_better_than_grid_best`)
- Topology collapse: `packages/treetime/src/commands/optimize/run.rs` (`find_zero_optimal_internal_edges`)
- Unimodal flag: `packages/treetime/src/gtr/gtr.rs` (`unimodal_branch_likelihood`)
- Related: [M-optimize-brent-zero-branch-non-unimodal](M-optimize-brent-zero-branch-non-unimodal.md) -- the boundary check resolves this defect
