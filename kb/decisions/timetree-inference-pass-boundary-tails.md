# Distribution Support-Boundary Tails and Inference-Pass Extrapolation

This document records an intentional deviation from v0 in how v1 handles evaluation of a discretized distribution outside its grid support, and how the timetree backward and forward passes assign per-side tail policies to the messages they produce.

## Background

A `GridFn` is a piecewise-linear function on a finite uniform grid. `Distribution::Function` reuses it as a finite-support probability distribution. These two roles impose conflicting requirements outside the grid: a generic interpolant has no defined value there, while a bounded probability density is zero there.

### v0 behavior

v0 stores neg-log probability and builds the interpolator with a constant fill value: `interp1d(xvals, yvals, kind='linear', fill_value=BIG_NUMBER, bounds_error=False)` at `packages/legacy/treetime/treetime/distribution.py:216`. `BIG_NUMBER` is `1e10` (`config.py:3`). Because the stored axis is neg-log probability, a fill of `1e10` corresponds to `exp(-1e10)`, effectively zero probability. So a v0 distribution is zero outside its support, and this is a property of the neg-log representation, not a generic extrapolation rule.

v0 does not express convolution-tail behavior as generic distribution extrapolation. `NodeInterpolator.convolve_fft()` at `packages/legacy/treetime/treetime/node_interpolator.py:231-256` separately constructs slope-based (linear in neg-log, exponential in probability) tails from the outermost trusted points, extending them only when the tail decays away from support and the domain boundary is not yet reached. There is no per-pass tail override.

### v1 problem

v1 `GridFn` previously did constant extrapolation unconditionally (return the nearest boundary value for any out-of-support query). Every distribution operation inherited that, so a bounded density silently gained an unbounded flat tail, diverging from v0 and from finite-support multiplication and division semantics.

## Decision

### Boundary policy on GridFn

`GridFn` carries an explicit per-side tail policy `BoundaryBehavior { Error, Zero, Constant }` in `left_extrap` and `right_extrap` (`packages/treetime-grid/src/grid_fn.rs`), defaulting to `Error`.

- `Error` (default): out-of-support evaluation returns an error. A bare grid function is a generic interpolant; a query outside its support is a programming error unless the caller has declared a tail. This makes out-of-support handling consistent across all `Distribution` variants (`Point` and `Range` already error out of support).
- `Zero`: return `0.0`. Under the plain-probability representation this is zero probability, matching the v0 result on the plain axis.
- `Constant`: return the nearest boundary value (flat tail). Used where a message is genuinely uninformative beyond an edge.

`GridFn` remains a policy-agnostic numeric primitive: `Zero` writes the literal `0.0`. That is only correct under a plain representation. Under the neg-log representation zero probability is `+inf`, not `0.0`, so a `Zero` tail on a neg-log `DistributionFunction` is rejected at construction (`DistributionFunction::with_left_extrap` / `with_right_extrap`, guarded by `YAxisPolicy::supports_zero_boundary`). The plain path is the one exercised by the timetree passes.

### Per-pass message tails

The backward and forward passes assign different tail policies to the messages they compute, applied immediately after each message is formed:

| Pass                      | Left tail (far past) | Right tail (far future) |
| ------------------------- | -------------------- | ----------------------- |
| Backward (leaves to root) | `Constant`           | `Zero`                  |
| Forward (root to leaves)  | `Zero`               | `Constant`              |

Backward message (`packages/treetime/src/timetree/inference/backward_pass.rs`): `parent_message = child_time_dist (*) (-branch_dist)`. The parent can be arbitrarily far in the past, so the left tail is `Constant`. The child's sampling date is a hard upper bound on the parent's age, so the right tail is `Zero`.

Forward message (`packages/treetime/src/timetree/inference/forward_pass.rs`): `dist_from_parent = parent_except_subtree (*) branch_dist`. The parent's time is a hard lower bound (non-negative branch lengths), so the left tail is `Zero`. There is no upper bound from the parent side on how far in the future the node could be, so the right tail is `Constant`.

`normalize()` rebuilds the grid and resets the tail policy to the default, so the forward pass applies the policy after normalization.

## Divergence from v0

This differs from v0 in two ways:

1. Tail ownership. v0 owns convolution tails inside `convolve_fft` as slope-based extrapolation; it has no per-pass tail override. v1 instead attaches a per-pass `Constant`/`Zero` policy to each message and does not reconstruct v0's slope-based convolution tails.
2. Default out-of-support handling. v0 returns effectively zero probability everywhere outside support (a soft value). v1 errors by default and requires an explicit `Zero` or `Constant` opt-in.

The v1 pass tails are motivated by node-time monotonicity: without a `Constant` left tail on the backward message, a parent's inferred time can be truncated too recent, producing negative branch lengths. The `Constant` left tail lets the combined distribution place the parent appropriately older. This is paired with, but distinct from, the child-time monotonicity clamp in the forward pass.

## Scope and interaction

Under v1's intersection-based multiplication, the product is evaluated only on the support intersection of its operands, so operand tails are not exercised by multiplication. Division also uses the exact support intersection while the divisor has the default `Error` boundary behavior. An explicit divisor tail extends the division domain on that side to the dividend boundary, allowing inference-message tails to take effect without changing finite-support division. Exact endpoint contact produces a point distribution, and positive-width products preserve the finest input grid spacing.

## Code references

- `packages/treetime-grid/src/grid_fn.rs`: `BoundaryBehavior`, `left_extrap`/`right_extrap`, `with_left_extrap`/`with_right_extrap`/`with_extrap`, `interp` (fallible), `resample` (propagates policy), `negate_arg_inplace` (swaps sides).
- `packages/treetime-distribution/src/policy.rs`: `YAxisPolicy::supports_zero_boundary` (`Plain` true, `NegLog` false).
- `packages/treetime-distribution/src/distribution_core/function.rs`: builder rejection of `Zero` under a non-supporting representation; `resample_dx` regrid boundary handling.
- `packages/treetime-distribution/src/distribution_core/distribution.rs`: `Distribution::with_left_extrap`/`with_right_extrap` (no-op for non-`Function` variants).
- `packages/treetime/src/timetree/inference/backward_pass.rs`, `forward_pass.rs`: per-pass tail application.
