# Distribution boundary tails and tail-aware arithmetic

**Type**: Intentional divergence from v0.

v1 discretizes probability distributions on finite uniform grids. Two questions arise: what value does a distribution take outside its grid, and how do arithmetic operations (multiplication, division) compute the result grid when operands carry explicit out-of-support policies. This document records both.

## Out-of-support evaluation

`GridFn` is a piecewise-linear interpolant on a finite uniform grid $[x_{\min}, x_{\max}]$. `Distribution::Function` reuses it as a probability density. These two roles disagree outside the grid: a generic interpolant has no defined value there, while a bounded probability density is zero. Convolution results in the timetree passes have a third behavior -- one side is unbounded but not zero (the parent could be arbitrarily older), so the density continues past the grid edge under a fitted decaying tail law.

### `BoundaryBehavior`

`GridFn` carries an independent tail policy on each side via `enum BoundaryBehavior` in `left_extrap` and `right_extrap`, defaulting to `Error`:

| Variant        | Class | Out-of-support value                                        | Use case                                                                                                                                                                                                                     |
| -------------- | ----- | ----------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Error`        | hard  | returns an error                                            | Default. A bare grid function is a generic interpolant; querying outside support is a programming error                                                                                                                      |
| `Hard`         | hard  | zero probability ($0$ under plain, $+\infty$ under neg-log) | The grid edge is the hard boundary: probability is zero beyond it, with no sub-grid gap to interpolate                                                                                                                       |
| `HardApproach` | hard  | fitted `HardApproachLaw` across the sub-grid gap            | The hard boundary lies below the first grid point; between them the density follows an edge-relative power-law-plus-linear approach law. Probability is still zero past the boundary itself                                  |
| `Constant`     | soft  | returns the boundary $y$                                    | Flat tail: genuinely uninformative beyond the edge. Non-integrable (infinite mass), so retained only for edges outside inference                                                                                             |
| `Linear`       | soft  | returns $p_{\text{edge}}\,e^{-k(t-t_{\text{edge}})}$        | Log-linear tail: a decaying exponential fitted from the edge points, a straight line in $-\ln p$. Carries the single slope $k$, has finite mass, and does not corrupt the quantile and HPD integrals the way `Constant` does |

Every variant is a complete value: the two law-carrying variants (`HardApproach`, `Linear`) always hold a fitted law, and the nullary variants (`Error`, `Hard`, `Constant`) declare a tail that needs none. A grid too small to fit a required law is an error at the fitting site, never a silent flat fallback, so no variant stands for "a law was wanted here but is missing".

A boundary is **hard** when the grid edge is a fact about the distribution -- probability is zero beyond (`Hard`, `HardApproach`), or evaluation beyond is undefined (`Error`). A boundary is **soft** when the grid edge is only where interpolation stopped and the distribution continues past it under a declared tail law (`Constant` or `Linear`). `fn BoundaryBehavior.is_soft()` is the predicate the arithmetic keys off: a soft boundary extends the evaluable domain, a hard boundary terminates it. Both soft tails route through this one predicate, so they share the arithmetic rules below.

`Linear` stores only the neg-log slope $k$ and re-reads the live grid edge on evaluation. A soft edge is a movable representation choice -- re-windowing and resampling shift it -- so anchoring the tail to the current edge keeps it valid across regridding, where a stored absolute anchor would go stale. The `Hard` approach law is edge-relative in the same way: only its boundary _location_ $t_\text{hard}$ is an immovable physical fact, while the ordinate is read from the live grid edge and the law stores only its shape (the power-law exponent and the linear slope), so it too survives regridding without a stored anchor. Both laws are shift-invariant: adding a constant to every $-\ln p$ leaves $k$, and the hard law's exponent and slope, unchanged.

`GridFn` is representation-agnostic: `Hard` writes the literal $0.0$, and the distribution layer maps that to the policy-correct zero-probability value ($0$ under plain, $+\infty$ under negative-log). Every variant is therefore valid under both axes. The timetree passes store distributions on the negative-log axis (`Distribution<NegLog>`), where the ordinate is $-\ln p$ and normalization is a pure ordinate shift.

Builder methods `fn GridFn.with_left_extrap()`, `fn GridFn.with_right_extrap()`, and `fn GridFn.with_extrap()` set the policy. `fn GridFn.resample()` propagates it. `fn GridFn.negate_arg_inplace()` swaps left and right (negating the argument reflects the domain).

## Tail assignments in the timetree pipeline

The inference pipeline produces several distribution types. Each carries per-side tails that downstream operations use when computing support intersections. Tails are metadata set on the convolution result; they are consumed by the next operation in the pipeline (multiplication or division).

| Distribution            | Left tail (past)  | Right tail (future) | Where set                                                                                                 | Consumed by                                                                     |
| ----------------------- | ----------------- | ------------------- | --------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| Leaf time constraint    | `Error`           | `Error`             | default                                                                                                   | backward convolution input                                                      |
| Branch length           | `HardApproach`    | `Linear` (fitted)   | `fn compute_branch_length_distribution()` (hard approach near `t=0`, `SoftTailLaw::fit` near `t_max`)     | convolution input (negated for backward); `fn compute_positional_log_lh()` eval |
| Backward parent message | `Linear` (fitted) | `Hard`              | `fn fit_message_soft_tail()` on the left, `with_right_extrap(Hard)` (backward_pass)                       | multiplication (combining children), division (forward cavity)                  |
| Internal node time dist | `Linear` (fitted) | `Hard`              | product result (left tail composed in closed form by `fn distribution_product()`, preserved by normalize) | forward convolution input, division dividend                                    |
| Forward message         | `Hard`            | `Linear` (fitted)   | `with_left_extrap(Hard)`, `fn fit_message_soft_tail()` on the right (forward_pass)                        | multiplication (refining node dist)                                             |
| Refined node time dist  | `Hard`            | `Hard`              | multiplication result (composed from forward message and subtree dist, preserved by normalize)            | `fn Distribution.likely_time()` extraction                                      |

### Why backward messages have `Linear` left / `Hard` right

The backward parent message is the convolution $\text{parent\_message} = \text{child\_time\_dist} \circledast (-\text{branch\_dist})$. It represents "when could the parent have existed, given this child?" The parent can be arbitrarily far in the past -- no child constrains how old its ancestor might be. The left side is therefore soft. It carries a fitted `Linear` tail rather than a flat `Constant`: the message genuinely decays away from its peak, and a log-linear tail captures that decay with finite mass, whereas a flat tail has infinite mass and corrupts the quantile and HPD integrals. The child's sampling date provides a hard upper bound on the parent's age (a parent cannot be more recent than its child), so the right tail is `Hard`.

### Why forward messages have `Hard` left / `Linear` right

The forward message is the convolution $\text{dist\_from\_parent} = \text{parent\_except\_subtree} \circledast \text{branch\_dist}$. It represents "when could this node have existed, given its parent and branch?" The parent's committed time provides a hard lower bound (branch lengths are non-negative), so the left tail is `Hard`. There is no upper bound from the parent side on how far in the future the node could be, so the right side is soft: a fitted `Linear` tail, for the same integrability reason as the backward left tail.

### Why the branch-length distribution has `HardApproach` left / `Linear` right

The branch-length distribution is the per-edge likelihood $L(t) \propto (\mu t)^{n} e^{-\mu t}$ over the branch duration, gridded on $[t_\text{min}, t_\text{max}]$. Its two edges are asymmetric facts. The left edge $t = 0$ is a hard boundary: a duration cannot be negative. For a branch with mutations the density vanishes there as a power law, so the left tail is a `HardApproach` law fitted near $t = 0$, which also carries the finite mode of a zero-mutation branch sitting on the boundary. The right edge $t_\text{max}$ is not a fact about the distribution, only where gridding stopped: the exponential $e^{-\mu t}$ keeps the density decaying past it (the Poisson indel term only adds a further linear neg-log slope), so the right side is a soft `Linear` tail fitted from the outermost points, with finite mass. Leaving it `Error` would declare a hard cutoff and drop any edge whose inferred duration exceeds $t_\text{max}$.

Convolution reads only the grid ordinates and reconstructs the message tails from its own output, so it does not consume these tails; `negate()` swaps and sign-flips them for the backward message. The tails are consumed directly by `fn compute_positional_log_lh()`, which evaluates the distribution at the inferred duration $\text{child\_time} - \text{parent\_time}$: a duration beyond $t_\text{max}$ now reads the extrapolated `Linear` tail instead of failing the `Error` default and dropping the edge.

### `normalize()` preserves tails

Under the negative-log axis, `fn Distribution.normalize()` subtracts the minimum ordinate via `fn DistributionFunction.shift_y()`. A shift is closed form on both fitted laws -- the soft-tail slope and the hard approach law's shape are shift-invariant -- so `shift_y` carries the laws through unchanged rather than dropping them. Every regrid also preserves the per-side policy and fitted law through the centralized `fn GridFn.regridded()` helper, and `fn GridFn.scale_y()` rescales the fitted laws rather than discarding them. Only an arbitrary `fn GridFn.mapv()` drops the law (keeping the side's class), which no normalization path uses.

This lets the backward pass combine child messages with a single product and normalization: the product composes the result tails (see [multiplication tails](#multiplication)) and normalization preserves them, so the folded result keeps its fitted `Linear` left tail and can still extend to reach a child with a disjoint finite grid. The child fold is `fn distribution_product()` (treetime-distribution), which co-locates the `Function` messages on one common grid and composes the summed left tail in closed form from the operands' own tails (`fn compose_multiplication_tail()`), so no tail is re-fit from the combined grid.

The forward pass applies its own tail policies _after_ normalization because the forward message needs `Hard` left and a fitted `Linear` right regardless of the tails the multiplication produced. It fits the right soft tail from the normalized grid; the shift left the slope unchanged, so fitting before or after normalization gives the same law.

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

| Operand tail | Behavior                                 | Rationale                                                                                                                                                |
| ------------ | ---------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Constant`   | Extend to the other operand's grid bound | The operand is evaluable at its boundary value beyond the grid; the product reflects the informative operand's shape scaled by a flat constant           |
| `Linear`     | Extend to the other operand's grid bound | The operand is evaluable at its exponential tail value beyond the grid; the product reflects the informative operand's shape scaled by the decaying tail |
| `Hard`       | Keep grid boundary                       | The value outside support is zero; $0 \cdot f(x) = 0$, so the product carries no information in the extension region                                     |
| `Error`      | Keep grid boundary                       | Out-of-support evaluation is undefined; strict finite-grid intersection                                                                                  |

Equivalently, per side the result takes the tightest (innermost) hard bound and the loosest (outermost) soft bound; a hard bound dominates a soft bound on the same side. Extending each soft side to the other operand's grid bound and intersecting selects this automatically.

The composed result tail per side is the strongest class of the two operands, `soft < hard < Error`, where the hard class covers both `Hard` and `HardApproach`. Two `Linear` tails compose in closed form: multiplication is addition in $-\ln p$, so their slopes add ($k = k_a + k_b$). A `Linear` tail times a flat `Constant` keeps the `Linear` slope, since a flat tail contributes slope zero. Two `HardApproach` laws compose their approach shapes likewise. `Linear` always carries a fitted slope; the total enum has no unfitted `Linear`, so no composition step defers a refit.

When both operands have `Error` or `Hard` tails (the default), the result is identical to strict intersection. When operands carry soft (`Constant` or `Linear`) tails -- as backward parent messages do (see [tail assignments](#tail-assignments-in-the-timetree-pipeline)) -- disjoint finite grids overlap via the tails. This prevents the product from collapsing to `Empty` when temporal signals conflict, for example under `--keep-root` where rerooting cannot resolve the tension between subtrees.

Point * Function follows the same per-side rule without constructing a result grid. A point in a `Constant` tail evaluates against the Function's boundary value, and a point in a `Linear` tail against the tail's exponential value. A point beyond an `Error` or `Hard` boundary produces `Empty`.

##### Empty invariant

A Function-producing multiplication returns `Empty` only when the operands' hard domains are genuinely disjoint. The hard domain of an operand is its grid bounds on hard sides and unbounded ($\pm\infty$) on soft sides, since a soft side continues under its tail law rather than terminating. Because a soft side always bridges a gap to the other operand, an empty product can arise only from two hard bounds facing each other across a gap; any other empty result is a numerical or logic collapse that would silently poison every ancestor to the root (the original motivating defect). `fn multiplication_empty_result()` checks the hard domains and reports `fn make_internal_error!()` rather than returning `Empty` when they overlap.

##### Result tails

A `Function` result carries per-side tails composed from the two operands' tails on that side. Beyond a boundary the product is evaluated pointwise: if either operand is undefined there the product is undefined (`Error`); otherwise if either operand is zero the product is zero (a hard class); only when both operands continue softly is the product soft. This is the maximum over the precedence `soft` < `hard` < `Error`, where soft covers `Constant`/`Linear` and hard covers `Hard`/`HardApproach` (the more restrictive tail wins):

| A tail  | B tail  | Result tail                                            |
| ------- | ------- | ------------------------------------------------------ |
| soft    | soft    | soft (`Linear` if either is `Linear`, else `Constant`) |
| soft    | hard    | hard                                                   |
| soft    | `Error` | `Error`                                                |
| hard    | hard    | hard                                                   |
| hard    | `Error` | `Error`                                                |
| `Error` | `Error` | `Error`                                                |

A `Range` or `Formula` operand has no interpolated tail (it is `Error` on both sides), so Range * Function and Formula * Function results carry `Error` tails. Only Function * Function can produce non-`Error` result tails, and only when both operands opt in. Combined with tail-preserving normalization, this is what lets the backward pass accumulate child messages without re-applying tails.

#### Division

The quotient spans the dividend's bounds intersected with the divisor's evaluable domain; the dividend's bounds are used as-is. This asymmetry reflects the use case: the forward pass divides a parent's time distribution (dividend) by a child's backward message (divisor) to compute the cavity. The divisor is never sampled beyond its real, grid-backed support: only a flat `Constant` tail extends the quotient past the divisor's grid edge, because it is the one tail evaluable there without distorting the quotient. Every other tail bounds the quotient at the divisor's grid edge, and the dividend's own tail carries the cavity beyond.

| Divisor tail            | Behavior                                          | Rationale                                                                                                                                                                                                                                                                                                      |
| ----------------------- | ------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Constant`              | Extend divisor bound to the dividend's grid bound | The divisor is evaluable at its flat boundary value, so dividing by it beyond the grid is well-defined and does not inflate the quotient; the cavity continues under the dividend                                                                                                                              |
| `Hard` / `HardApproach` | Bound the quotient at the divisor's grid edge     | Probability is zero beyond the boundary. Under neg-log that is $+\infty$, so `dividend - (+\infty) = -\infty` -- a spurious spike that collapses the downstream convolution. `fn YAxisPolicy::safe_divisor()` cannot rescue it: it floors small divisors under plain storage but is the identity under neg-log |
| `Linear`                | Bound the quotient at the divisor's grid edge     | The divisor decays beyond the grid, so dividing by its extrapolated tail inflates the quotient into a spurious spike                                                                                                                                                                                           |
| `Error`                 | Bound the quotient at the divisor's grid edge     | Out-of-support evaluation is undefined; finite-support division constrained to the divisor's grid                                                                                                                                                                                                              |

This is why an independently re-gridded dividend (for example, a mass-sized parent posterior whose grid extends past a hard-bounded child message) no longer collapses the cavity: the quotient stops at the divisor's grid edge instead of sampling the divisor's out-of-support extrapolation.

Coalescent contributions do not introduce another grid. The backward pass multiplies child messages first and evaluates the coalescent contribution pointwise on the resulting grid in negative-log space.

## Divergence from v0

Two specific divergences from v0's distribution handling:

1. **Convolution tail reconstruction.** v0 rebuilds convolution tails inside `fn NodeInterpolator.convolve_fft()` ([node_interpolator.py#L231-L256](../../packages/legacy/treetime/treetime/node_interpolator.py#L231-L256)) as slope-based extrapolation (linear in neg-log, exponential in probability) from the outermost trusted points. v1 now reconstructs the same tails in `fn convolution_function_function()`: the FFT runs in plain probability space and is trusted only above a peak-relative floor (`1e-13`, matching v0), and the sub-floor tail is rebuilt by a two-point secant slope in negative-log space, extended only where it decays away from support. The reconstruction is policy-generic (`Plain` and `NegLog`) because each operand is converted to peak-normalized plain probability around the FFT via `fn YAxisPolicy::to_neg_log()` / `fn YAxisPolicy::from_neg_log()` and back. The inference passes then fit a per-message soft tail from the reconstructed grid (`fn fit_message_soft_tail()`) and attach it as the message's out-of-support policy on the soft side, keeping the hard side `Hard`. That policy governs how the _following_ multiplication or division extends the support intersection.

2. **Default out-of-support handling.** v0 returns effectively zero probability outside support as a soft value (`fill_value=1e10` in neg-log, i.e. $\exp(-10^{10})$). v1 returns `Error` by default, requiring an explicit tail opt-in.

The v1 pass soft tails serve two ends. They preserve node-time monotonicity -- without a soft left tail on the backward message a parent's inferred time can be truncated too recent, producing negative branch lengths, and the soft tail lets the combined distribution place the parent appropriately older. The fitted `Linear` tail also has finite mass, so it does not corrupt the quantile and HPD integrals the way the earlier flat `Constant` tail did; because log-linear extrapolation slightly over-estimates the far tail it can never manufacture a spurious zero, so the monotonicity guarantee is retained. A separate forward-pass projection currently clamps committed internal-node point estimates; the statistical contract for that projection remains open in [kb/issues/M-timetree-marginal-node-times-can-violate-topology.md](../issues/M-timetree-marginal-node-times-can-violate-topology.md).

The hard side of each message currently uses the nullary `Hard` variant, which places zero probability past the grid edge. Fitting a `HardApproach` law there (to carry a mode sitting on the hard boundary, as for a zero-mutation branch) is a separate, unresolved modeling decision tracked in [kb/issues/N-grid-hard-approach-law-drops-linear-and-folds-indels.md](../issues/N-grid-hard-approach-law-drops-linear-and-folds-indels.md).

## Accepted limitations

This is a bounded resampling rule, not a derived interpolation-error criterion. Rounding can make the realized spacing slightly finer or coarser than the target. The $10^6$-point cap can make it substantially coarser. Sequential multiplication resamples at each binary operation and is not guaranteed to be associative, although Function * Function multiplication is commutative because its support and spacing choices are symmetric.

Convolution now routes both policies through a shared plain-space peak-normalization (`fn to_peak_normalized_plain()`), so the `Plain` result carries ULP-level round-trip noise from the `-ln`/`exp` bracketing that the previous raw-probability path did not. The difference is at the limit of `f64` precision, far below any scientific tolerance, but exact-equality assertions on `Plain` convolution output no longer hold.

Scientific workflows requiring posterior peak, normalization, or integrated-probability error bounds need a separate accuracy contract.

## Code references

### Tail policy

- [packages/treetime-grid/src/boundary_behavior.rs](../../packages/treetime-grid/src/boundary_behavior.rs): `enum BoundaryBehavior` (`Error`/`Hard`/`HardApproach`/`Constant`/`Linear`, every variant a complete value), `fn BoundaryBehavior.is_soft()`/`is_hard()`, `const DEFAULT_TAIL_FIT_POINTS`
- [packages/treetime-grid/src/soft_tail_law.rs](../../packages/treetime-grid/src/soft_tail_law.rs): `struct SoftTailLaw` (single neg-log slope), `fn SoftTailLaw::fit()`, `fn SoftTailLaw::mass()`, `fn SoftTailLaw::compose_multiply()`
- [packages/treetime-grid/src/grid_fn.rs](../../packages/treetime-grid/src/grid_fn.rs): `left_extrap`/`right_extrap` fields, `fn GridFn.with_left_extrap()`/`fn GridFn.with_right_extrap()`/`fn GridFn.with_extrap()`, `fn GridFn.interp()` (fallible, dispatches by tail), `fn GridFn.regridded()` (centralized policy/law carry across every regrid), `fn GridFn.resample()`, `fn GridFn.scale_y()` (rescales fitted laws), `fn GridFn.shift_y()` (shift-invariant, carries laws), `fn GridFn.mapv()` (drops the law, keeps the class), `fn GridFn.negate_arg_inplace()` (swaps sides)
- [packages/treetime-distribution/src/distribution_core/distribution.rs](../../packages/treetime-distribution/src/distribution_core/distribution.rs): `fn Distribution.with_left_extrap()`/`fn Distribution.with_right_extrap()` (no-op for non-Function variants); `fn Distribution.normalize()` (subtracts the minimum ordinate via `fn DistributionFunction.shift_y()` under `NegLog`)

### Tail-aware arithmetic

- [packages/treetime-distribution/src/distribution_ops/multiply.rs](../../packages/treetime-distribution/src/distribution_ops/multiply.rs): `fn multiplication_support_intersection()`, `fn multiplication_empty_result()` (empty invariant), `fn hard_domains_disjoint()`, `fn compose_multiplication_tail()`, `fn with_composed_tails()`, `fn multiply_point_function()`, `fn multiply_function_function()`, `fn multiply_range_function()`
- [packages/treetime-distribution/src/distribution_ops/divide.rs](../../packages/treetime-distribution/src/distribution_ops/divide.rs): `fn division_support_intersection()`, `fn divide_function_by_function()`, `fn divide_range_by_function()`

### Inference pass tail application

- [packages/treetime/src/timetree/inference/tail_fit.rs](../../packages/treetime/src/timetree/inference/tail_fit.rs): `fn fit_message_soft_tail()` -- fits the soft-side `Linear` tail of a message, errors on a degenerate grid, inert on non-Function messages. Shared by both passes.
- [packages/treetime/src/timetree/inference/backward_pass.rs](../../packages/treetime/src/timetree/inference/backward_pass.rs): backward message soft-left / hard-right assignment; the child fold delegates to `fn distribution_product()` (treetime-distribution), which composes the summed tails in closed form.
- [packages/treetime-distribution/src/distribution_ops/multiply.rs](../../packages/treetime-distribution/src/distribution_ops/multiply.rs): `fn distribution_product()` folds the child `Function` messages on one common grid; `fn compose_multiplication_tail()` composes each side's result tail in closed form (soft slopes add, a hard bound dominates).
- [packages/treetime/src/timetree/inference/forward_pass.rs](../../packages/treetime/src/timetree/inference/forward_pass.rs): forward message hard-left / soft-right assignment, fit after `normalize()`.

### Tests

- [packages/treetime-distribution/src/distribution_ops/**tests**/test_multiply.rs](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_multiply.rs): tail combinations (overlapping, disjoint, mixed `Constant`/`Hard`/`Error`, commutativity)
- [packages/treetime-distribution/src/distribution_ops/**tests**/test_divide.rs](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_divide.rs): exact intersection under default boundaries, extension under explicit divisor tails
- [packages/treetime-distribution/src/distribution_core/**tests**/test_boundary_behavior.rs](../../packages/treetime-distribution/src/distribution_core/__tests__/test_boundary_behavior.rs): boundary-variant evaluation and composition
- [packages/treetime/src/timetree/inference/**tests**/test_tail_fit.rs](../../packages/treetime/src/timetree/inference/__tests__/test_tail_fit.rs): `fn fit_message_soft_tail()` -- slope recovery, finite tail mass, inert on non-Function, error on degenerate grid

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
