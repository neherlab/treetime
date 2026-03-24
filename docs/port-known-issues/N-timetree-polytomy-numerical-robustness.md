# Polytomy resolution numerical robustness

Several defensive programming gaps in polytomy resolution code. None are active bugs in the current call sequence, but they create latent defects that would activate if the call sequence changes or inputs are degenerate.

v1: [`packages/treetime/src/commands/timetree/optimization/polytomy.rs`](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs)

## Items

### `ln(0.0)` produces `-Inf` in cost function

`MergeCostFunction::cost()` computes `.unwrap_or(1e-10).ln()` on `Distribution::eval()` results. The `unwrap_or` guard catches `Err` but not `Ok(0.0)`. `Distribution::Function` variants can return `Ok(0.0)` when boundary values underflow to zero during `exp()` normalization in `compute_branch_length_distribution()`. The proper fix is clamping `normalized_prob` at the source in [`packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L53`](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L53), not a consumer-side clamp.

v0 avoids this by working in neg-log space with `fill_value=BIG_NUMBER`.

### `Distribution::Empty` fallback silently prevents merges

Missing `branch_length_distribution` defaults to `Distribution::Empty`, which always errors on `eval()`. Cost gain reduces to `-penalty` (always negative), silently preventing resolution. The current call sequence guarantees distributions are populated before polytomy resolution, but a `debug_assert!` or log warning would catch future ordering changes.

### Undocumented `1e-10` fallback probability

The `unwrap_or(1e-10)` value is an arbitrary floor with no derivation from model parameters. When both branches hit the fallback, the effect cancels. In the mixed case (one branch valid, one fallback), the floor contributes `ln(1e-10) = -23` to cost.

### Fixed `1e-10` Brent bounds offset ignores time scale

`BrentOpt::new(parent_time + 1e-10, child_min_time - 1e-10)` uses a fixed offset. At time values around `1e6` (deep phylogenetics), machine epsilon is `~1.2e-10`, making the offset barely representable. v0 scales tolerance via `1e-4 * one_mutation`. TreeTime targets viral phylogenetics (times 1900-2025) where `1e-10` is safely representable.

### No postcondition on branch length non-negativity

`merge_children()` computes new edge `time_length` values without asserting non-negativity. The Brent optimizer bounds guarantee positive branch lengths under normal conditions, but floating-point edge cases at extreme time scales are unguarded. A `debug_assert!(length >= 0.0)` would catch regressions.

### Silent fallback on optimizer failure

`compute_merge_gain()` catches Brent optimizer failures and silently falls back to midpoint evaluation. The fallback maps errors to zero gain (preventing merge), which is benign. Adding `debug!` logging would aid diagnosis.

### `unwrap_or(0.0)` on missing node times

`collect_children_info()` and `resolve_single_polytomy()` use `time.unwrap_or(0.0)` for `None` node times. In calendar time, 0.0 is year 0 CE, producing extreme branch lengths. The current call sequence guarantees all nodes have times before polytomy resolution, but the fallback masks ordering bugs.
