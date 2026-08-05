# Distribution boundary tails and tail-aware arithmetic

**Type**: Intentional divergence from v0.

v1 discretizes probability distributions on finite uniform grids. Two questions arise: what value does a distribution take outside its grid, and how do arithmetic operations (multiplication, division) compute the result grid when operands carry explicit out-of-support policies. This document records both.

## Out-of-support evaluation

`GridFn` is a piecewise-linear interpolant on a finite uniform grid $[x_{\min}, x_{\max}]$. `Distribution::Function` reuses it as a probability density. These two roles disagree outside the grid: a generic interpolant has no defined value there, while a bounded probability density is zero. Convolution results in the timetree passes have a third behavior -- one side may be genuinely uninformative (the parent could be arbitrarily older), so the distribution is flat, not zero, beyond the grid edge.

### `BoundaryBehavior`

`GridFn` carries an independent tail policy on each side via `enum BoundaryBehavior` in `left_extrap` and `right_extrap`, defaulting to `Error`:

| Variant    | Out-of-support value     | Use case                                                                                                                |
| ---------- | ------------------------ | ----------------------------------------------------------------------------------------------------------------------- |
| `Error`    | returns an error         | Default. A bare grid function is a generic interpolant; querying outside support is a programming error                 |
| `Zero`     | returns $0.0$            | Bounded probability density: zero probability outside support. Matches v0's result on the plain-probability axis        |
| `Constant` | returns the boundary $y$ | Flat tail: the distribution is genuinely uninformative beyond the grid edge, so the boundary value extends indefinitely |

`GridFn` is representation-agnostic: `Zero` writes the literal $0.0$. Under a neg-log representation, zero probability is $+\infty$, not $0.0$, so `fn DistributionFunction.with_left_extrap()` and `fn DistributionFunction.with_right_extrap()` reject `Zero` when the representation does not support it (guarded by `fn YAxisPolicy::supports_zero_boundary()`). The plain-probability path used by the timetree passes accepts all three variants.

Builder methods `fn GridFn.with_left_extrap()`, `fn GridFn.with_right_extrap()`, and `fn GridFn.with_extrap()` set the policy. `fn GridFn.resample()` propagates it. `fn GridFn.negate_arg_inplace()` swaps left and right (negating the argument reflects the domain).

## Tail assignments in the timetree pipeline

The inference pipeline produces several distribution types. Each carries per-side tails that downstream operations use when computing support intersections. Tails are metadata set on the convolution result; they are consumed by the next operation in the pipeline (multiplication or division).

| Distribution            | Left tail (past) | Right tail (future) | Where set                                                                                               | Consumed by                                                    |
| ----------------------- | ---------------- | ------------------- | ------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------- |
| Leaf time constraint    | `Error`          | `Error`             | default                                                                                                 | backward convolution input                                     |
| Branch length           | `Error`          | `Error`             | default                                                                                                 | convolution input (negated for backward)                       |
| Backward parent message | `Constant`       | `Zero`              | [backward_pass.rs#L125-L126](../../packages/treetime/src/timetree/inference/backward_pass.rs#L125-L126) | multiplication (combining children), division (forward cavity) |
| Internal node time dist | `Constant`       | `Zero`              | multiplication result (composed from child messages, preserved by normalize)                            | forward convolution input, division dividend                   |
| Forward message         | `Zero`           | `Constant`          | [forward_pass.rs#L126-L127](../../packages/treetime/src/timetree/inference/forward_pass.rs#L126-L127)   | multiplication (refining node dist)                            |
| Refined node time dist  | `Zero`           | `Zero`              | multiplication result (composed from forward message and subtree dist, preserved by normalize)          | `fn Distribution.likely_time()` extraction                     |

### Why backward messages have `Constant` left / `Zero` right

The backward parent message is the convolution $\text{parent\_message} = \text{child\_time\_dist} \circledast (-\text{branch\_dist})$. It represents "when could the parent have existed, given this child?" The parent can be arbitrarily far in the past -- no child constrains how old its ancestor might be. The left tail is therefore `Constant`: the message is flat (uninformative) to the left of its grid. The child's sampling date provides a hard upper bound on the parent's age (a parent cannot be more recent than its child), so the right tail is `Zero`.

### Why forward messages have `Zero` left / `Constant` right

The forward message is the convolution $\text{dist\_from\_parent} = \text{parent\_except\_subtree} \circledast \text{branch\_dist}$. It represents "when could this node have existed, given its parent and branch?" The parent's committed time provides a hard lower bound (branch lengths are non-negative), so the left tail is `Zero`. There is no upper bound from the parent side on how far in the future the node could be, so the right tail is `Constant`.

### `normalize()` preserves tails

`fn Distribution.normalize()` calls `fn Distribution.scale_by()`, which calls `fn DistributionFunction.scale_y()`. `fn DistributionFunction.scale_y()` scales the y values with `fn GridFn.mapv()`, which copies both tail policies. Scaling by a positive factor does not change out-of-support behavior, so the tails survive normalization.

This lets the backward pass combine child messages with multiplication and normalization alone: multiplication composes the result tails (see [multiplication tails](#multiplication)) and normalization preserves them, so a chain of multiplications keeps its `Constant` left tail and can still extend to reach a child with a disjoint finite grid. No manual re-application is needed.

The forward pass still applies its own tail policies _after_ normalization ([forward_pass.rs#L126-L127](../../packages/treetime/src/timetree/inference/forward_pass.rs#L126-L127), [forward_pass.rs#L141-L142](../../packages/treetime/src/timetree/inference/forward_pass.rs#L141-L142)) because the forward message needs `Zero` left / `Constant` right regardless of the tails the multiplication produced. This explicit application overwrites the preserved tails, which is correct.

## Grid intersection contract

v1 computes the result grid for both multiplication and division by resampling onto a new uniform grid spanning the exact analytical support intersection. This is an intentional divergence from v0, which collects knots from both operands inside the overlap and builds a non-uniform grid from them ([distribution.py#L82-L145](../../packages/legacy/treetime/treetime/distribution.py#L82-L145)).

Given an intersection $[x_{\min}, x_{\max}]$, the result grid has:

$$n = \text{clamp}\left(\text{round}\left(\frac{x_{\max} - x_{\min}}{\Delta x}\right) + 1,\; 2,\; 10^6\right)$$

where $\Delta x = \min(\Delta x_a, \Delta x_b)$ for Function * Function, and $\Delta x = \Delta x_f$ for Range * Function or Formula * Function (using the Function operand's spacing).

The result grid includes both analytical endpoints via `Array1::linspace`. A disjoint intersection produces `Empty`. Endpoint-only contact produces a `Point` distribution evaluated at the shared boundary, matching v0's behavior of converting one surviving knot to a delta.

Endpoint contact uses exact comparison of the analytical bounds, not a tolerance. A tolerance would enlarge the intersection and could turn two narrowly disjoint supports into a `Point`.

The $10^6$-point cap prevents unbounded allocations when a narrow input spacing combines with a wide intersection.

### How tails extend the intersection

Both multiplication and division compute the intersection after extending each operand's evaluable domain based on its tail policy. The extension is per-side and independent: each side of each operand is resolved from its own `BoundaryBehavior`, and the result grid spans the intersection of the two extended domains. Without extension, operands with disjoint finite grids but overlapping via tails would collapse to `Empty`.

#### Multiplication

Each operand extends its own domain symmetrically. The extension limit is the other operand's grid bound on the same side:

| Operand tail | Behavior                                 | Rationale                                                                                                                                      |
| ------------ | ---------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `Constant`   | Extend to the other operand's grid bound | The operand is evaluable at its boundary value beyond the grid; the product reflects the informative operand's shape scaled by a flat constant |
| `Zero`       | Keep grid boundary                       | The value outside support is zero; $0 \cdot f(x) = 0$, so the product carries no information in the extension region                           |
| `Error`      | Keep grid boundary                       | Out-of-support evaluation is undefined; strict finite-grid intersection                                                                        |

When both operands have `Error` or `Zero` tails (the default), the result is identical to strict intersection. When operands carry `Constant` tails -- as backward parent messages do (see [tail assignments](#tail-assignments-in-the-timetree-pipeline)) -- disjoint finite grids overlap via the tails. This prevents the product from collapsing to `Empty` when temporal signals conflict, for example under `--keep-root` where rerooting cannot resolve the tension between subtrees.

Point * Function follows the same per-side rule without constructing a result grid. A point in a `Constant` tail evaluates against the Function's boundary value. A point beyond an `Error` or `Zero` boundary produces `Empty`.

##### Result tails

A `Function` result carries per-side tails composed from the two operands' tails on that side. Beyond a boundary the product is evaluated pointwise: if either operand is undefined there the product is undefined (`Error`); otherwise if either operand is zero the product is zero (`Zero`); only when both operands are flat non-zero constants is the product a flat constant (`Constant`). This is the maximum over the precedence `Constant` < `Zero` < `Error` (the more restrictive tail wins):

| A tail     | B tail     | Result tail |
| ---------- | ---------- | ----------- |
| `Constant` | `Constant` | `Constant`  |
| `Constant` | `Zero`     | `Zero`      |
| `Constant` | `Error`    | `Error`     |
| `Zero`     | `Zero`     | `Zero`      |
| `Zero`     | `Error`    | `Error`     |
| `Error`    | `Error`    | `Error`     |

A `Range` or `Formula` operand has no interpolated tail (it is `Error` on both sides), so Range * Function and Formula * Function results carry `Error` tails. Only Function * Function can produce non-`Error` result tails, and only when both operands opt in. Combined with tail-preserving normalization, this is what lets the backward pass accumulate child messages without re-applying tails.

#### Division

The divisor extends its domain; the dividend's bounds are used as-is. This asymmetry reflects the use case: the forward pass divides a parent's time distribution (dividend) by a child's backward message (divisor) to compute the cavity. The divisor's tail metadata determines how far the cavity extends:

| Divisor tail | Behavior                                          | Rationale                                                                               |
| ------------ | ------------------------------------------------- | --------------------------------------------------------------------------------------- |
| `Constant`   | Extend divisor bound to the dividend's grid bound | The divisor is evaluable at its boundary value; the cavity reflects the dividend        |
| `Zero`       | Extend divisor bound to the dividend's grid bound | The divisor is evaluable at zero; `fn YAxisPolicy::safe_divisor()` handles the division |
| `Error`      | Keep divisor grid boundary                        | Finite-support division, constrained to the divisor's grid                              |

Coalescent contributions do not introduce another grid. The backward pass multiplies child messages first and evaluates the coalescent contribution pointwise on the resulting grid in negative-log space.

## Divergence from v0

Two specific divergences from v0's distribution handling:

1. **Convolution tail reconstruction.** v0 rebuilds convolution tails inside `fn NodeInterpolator.convolve_fft()` ([node_interpolator.py#L231-L256](../../packages/legacy/treetime/treetime/node_interpolator.py#L231-L256)) as slope-based extrapolation (linear in neg-log, exponential in probability) from the outermost trusted points. v1 now reconstructs the same tails in `fn convolution_function_function()`: the FFT runs in plain probability space and is trusted only above a peak-relative floor (`1e-13`, matching v0), and the sub-floor tail is rebuilt by a two-point secant slope in negative-log space, extended only where it decays away from support. The reconstruction is policy-generic (`Plain` and `NegLog`) because each operand is converted to peak-normalized plain probability around the FFT via `fn YAxisPolicy::to_neg_log()` / `fn YAxisPolicy::from_neg_log()` and back. Separately, the inference passes still attach a per-message `Constant`/`Zero` out-of-support policy; that policy governs how the _following_ multiplication or division extends the support intersection and is independent of the on-grid tail values the convolution reconstructs.

2. **Default out-of-support handling.** v0 returns effectively zero probability outside support as a soft value (`fill_value=1e10` in neg-log, i.e. $\exp(-10^{10})$). v1 returns `Error` by default, requiring an explicit `Zero` or `Constant` opt-in.

The v1 pass tails are motivated by node-time monotonicity: without a `Constant` left tail on the backward message, a parent's inferred time can be truncated too recent, producing negative branch lengths. The `Constant` left tail lets the combined distribution place the parent appropriately older. A separate forward-pass projection currently clamps committed internal-node point estimates; the statistical contract for that projection remains open in [kb/issues/M-timetree-marginal-node-times-can-violate-topology.md](../issues/M-timetree-marginal-node-times-can-violate-topology.md).

## Accepted limitations

This is a bounded resampling rule, not a derived interpolation-error criterion. Rounding can make the realized spacing slightly finer or coarser than the target. The $10^6$-point cap can make it substantially coarser. Sequential multiplication resamples at each binary operation and is not guaranteed to be associative, although Function * Function multiplication is commutative because its support and spacing choices are symmetric.

Convolution now routes both policies through a shared plain-space peak-normalization (`fn to_peak_normalized_plain()`), so the `Plain` result carries ULP-level round-trip noise from the `-ln`/`exp` bracketing that the previous raw-probability path did not. The difference is at the limit of `f64` precision, far below any scientific tolerance, but exact-equality assertions on `Plain` convolution output no longer hold.

Scientific workflows requiring posterior peak, normalization, or integrated-probability error bounds need a separate accuracy contract.

## Code references

### Tail policy

- [packages/treetime-grid/src/grid_fn.rs](../../packages/treetime-grid/src/grid_fn.rs): `enum BoundaryBehavior`, `left_extrap`/`right_extrap` fields, `fn GridFn.with_left_extrap()`/`fn GridFn.with_right_extrap()`/`fn GridFn.with_extrap()`, `fn GridFn.interp()` (fallible, dispatches by tail), `fn GridFn.resample()` (propagates policy), `fn GridFn.negate_arg_inplace()` (swaps sides)
- [packages/treetime-distribution/src/policy.rs](../../packages/treetime-distribution/src/policy.rs): `fn YAxisPolicy::supports_zero_boundary()` (`Plain` true, `NegLog` false)
- [packages/treetime-distribution/src/distribution_core/function.rs](../../packages/treetime-distribution/src/distribution_core/function.rs): builder rejection of `Zero` under neg-log; `fn DistributionFunction.scale_y()` (preserves tails via `fn GridFn.mapv()`); `fn DistributionFunction.resample_dx()` regrid boundary handling
- [packages/treetime-distribution/src/distribution_core/distribution.rs](../../packages/treetime-distribution/src/distribution_core/distribution.rs): `fn Distribution.with_left_extrap()`/`fn Distribution.with_right_extrap()` (no-op for non-Function variants); `fn Distribution.normalize()` (preserves tails via `fn Distribution.scale_by()`)

### Tail-aware arithmetic

- [packages/treetime-distribution/src/distribution_ops/multiply.rs](../../packages/treetime-distribution/src/distribution_ops/multiply.rs): `fn multiplication_support_intersection()`, `fn compose_multiplication_tail()`, `fn with_composed_tails()`, `fn multiply_point_function()`, `fn multiply_function_function()`, `fn multiply_range_function()`
- [packages/treetime-distribution/src/distribution_ops/divide.rs](../../packages/treetime-distribution/src/distribution_ops/divide.rs): `fn division_support_intersection()`, `fn divide_function_by_function()`, `fn divide_range_by_function()`

### Inference pass tail application

- [packages/treetime/src/timetree/inference/backward_pass.rs](../../packages/treetime/src/timetree/inference/backward_pass.rs): backward message tail assignment (lines 125-126). Child combination is multiply + normalize only; tail composition and preservation remove the former manual re-application.
- [packages/treetime/src/timetree/inference/forward_pass.rs](../../packages/treetime/src/timetree/inference/forward_pass.rs): forward message tail assignment (lines 126-127, 141-142)

### Tests

- [packages/treetime-distribution/src/distribution_ops/**tests**/test_multiply.rs](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_multiply.rs): tail combinations (overlapping, disjoint, mixed `Constant`/`Zero`/`Error`, commutativity)
- [packages/treetime-distribution/src/distribution_ops/**tests**/test_divide.rs](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_divide.rs): exact intersection under default boundaries, extension under explicit divisor tails
- [packages/treetime-distribution/src/distribution_core/**tests**/test_boundary_behavior.rs](../../packages/treetime-distribution/src/distribution_core/__tests__/test_boundary_behavior.rs): representation constraints (`Zero` accepted by plain, rejected by neg-log; `Constant` accepted by both)

## Related knowledge base items

### Specifications

- [kb/algo/timetree-extrapolation-and-time-clamping.md](../algo/timetree-extrapolation-and-time-clamping.md): original spec defining extrapolation policies, exact intersection boundaries, and the monotonicity clamp. Source material for this decision.

### Proposals

- [kb/proposals/exponential-branch-length-tails.md](../proposals/exponential-branch-length-tails.md): a `BoundaryBehavior::Exponential` variant for branch length distributions, continuing the log-linear slope at the grid boundary. Would replace the flat `Constant` tail on backward messages with a physically motivated exponential decay.

### Issues

- [kb/issues/M-timetree-silent-empty-time-distribution.md](../issues/M-timetree-silent-empty-time-distribution.md): `fn set_likely_time()` silently skips nodes whose time distribution is `Empty`. Safety net for cases where tail-aware multiplication still cannot prevent empty results.
- [kb/issues/M-timetree-marginal-node-times-can-violate-topology.md](../issues/M-timetree-marginal-node-times-can-violate-topology.md): independent marginal modes can invert parent-child ordering. The forward-pass monotonicity clamp addresses this but its statistical contract is unresolved.
- [kb/issues/M-timetree-gm-runner-missing-internal-times.md](../issues/M-timetree-gm-runner-missing-internal-times.md): golden master tests for 5 datasets fail because internal node times are absent. Root cause (empty distributions from disjoint multiplication) is addressed by tail-aware multiplication; verification blocked by separate grid-width tolerance ignores.
- [kb/issues/M-distribution-normalization-erases-errors.md](../issues/M-distribution-normalization-erases-errors.md): `fn Distribution.normalize()` returns `Empty` for non-positive or non-finite max values, silently discarding error information.
- [kb/issues/N-distribution-function-product-negative-roundoff.md](../issues/N-distribution-function-product-negative-roundoff.md): multiplication of near-zero values can produce negative roundoff, which downstream consumers may not expect.
