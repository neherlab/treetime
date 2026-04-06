# Brent optimizer failures silently return stale branch length

## Problem

`brent_inner()` discards `argmin::BrentOpt` solver errors and returns the incoming `branch_length` unchanged:

```rust
match result {
  Ok(res) => res.state().best_param.unwrap_or(branch_length),
  Err(_) => branch_length,
}
```

A failed per-edge optimization is indistinguishable from a successful no-op. The caller (`run_optimize_mixed`) has no way to detect that the edge was not optimized.

## Impact

Bracket mistakes, evaluation-domain errors, or solver failures leave edges at their previous length while the caller believes Brent optimization succeeded. This biases inferred branch lengths and makes downstream likelihood comparisons unreliable. The failure is silent -- no log message, no error propagation.

## Approach

Two options:

1. Return `Result<f64, Report>` from `brent_inner()` and propagate solver failure through `run_optimize_mixed()`.
2. If a non-fatal fallback is needed: on `Err`, fall back to `grid_search_inner()` and emit a warning identifying the affected edge. This preserves forward progress while keeping the optimization path auditable.

The `Ok` path also silently falls back via `unwrap_or(branch_length)` when `best_param` is `None`. This should be treated as a failure too.

## Cross-references

- Brent error handling: `packages/treetime/src/commands/optimize/method_brent.rs` (`brent_inner`, match block)
- Caller: `packages/treetime/src/commands/optimize/optimize_unified.rs` (`run_optimize_mixed`)
