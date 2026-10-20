# Distribution boundary tails and tail-aware arithmetic

**Type**: Intentional divergence from v0.

v1 discretizes probability distributions on finite uniform grids. Two questions arise: what value does a distribution take outside its grid, and how do arithmetic operations (multiplication, division) compute the result grid when operands carry explicit out-of-support policies. This document records both.

## Out-of-support evaluation

`GridFn` is a piecewise-linear interpolant on a finite uniform grid $[x_{\min}, x_{\max}]$. `Distribution::Function` reuses it as a probability density. These two roles disagree outside the grid: a generic interpolant has no defined value there, while a bounded probability density is zero. Convolution results in the timetree passes have a third behavior -- one side is unbounded but not zero (the parent could be arbitrarily older), so the density continues past the grid edge under a fitted decaying tail law.

### `BoundaryBehavior`

`DistributionFunction` carries an independent tail policy on each side via `enum BoundaryBehavior` in `left_extrap` and `right_extrap`, defaulting to `Error`. `GridFn` remains a representation-agnostic interpolant; callers supply policies when they request extrapolation.

| Variant        | Class | Out-of-support value                                        | Use case                                                                                                                                                                                                                                                                         |
| -------------- | ----- | ----------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Error`        | hard  | returns an error                                            | Default for a Function distribution; querying outside support is a programming error                                                                                                                                                                                             |
| `Hard`         | hard  | zero probability ($0$ under plain, $+\infty$ under neg-log) | The grid edge is the hard boundary: probability is zero beyond it, with no sub-grid gap to interpolate                                                                                                                                                                           |
| `HardApproach` | hard  | fitted `HardApproachLaw` across the sub-grid gap            | The hard boundary lies below the first grid point; between them the density follows an edge-relative single-exponent power-law approach law. A finite mode on the boundary is instead stored as an exact grid endpoint (`Hard`), so this variant carries only the divergent case |
| `Linear`       | soft  | returns $p_{\text{edge}}\,e^{-k(t-t_{\text{edge}})}$        | Log-linear tail: a decaying exponential fitted from the edge points, a straight line in $-\ln p$. Carries the single slope $k$ and has finite mass, so the quantile and interval integrals stay finite                                                                                          |

Every variant is a complete value: the two law-carrying variants (`HardApproach`, `Linear`) always hold a fitted law, and the nullary variants (`Error`, `Hard`) declare a tail that needs none. A grid too small to fit a required law is an error at the fitting site, never a silent flat fallback, so no variant stands for "a law was wanted here but is missing".

A boundary is **hard** when the grid edge is a fact about the distribution -- probability is zero beyond (`Hard`, `HardApproach`), or evaluation beyond is undefined (`Error`). A boundary is **soft** when the grid edge is only where interpolation stopped and the distribution continues past it under its declared `Linear` tail law. `fn BoundaryBehavior.is_soft()` is the predicate the arithmetic keys off: a soft boundary extends the evaluable domain, a hard boundary terminates it.

`Linear` stores only the neg-log slope $k$ and re-reads the live grid edge on evaluation. A soft edge is a movable representation choice -- re-windowing and resampling shift it -- so anchoring the tail to the current edge keeps it valid across regridding, where a stored absolute anchor would go stale. The `Hard` approach law is edge-relative in the same way: only its boundary _location_ $t_\text{hard}$ is an immovable physical fact, while the ordinate is read from the live grid edge and the law stores only its shape (the single power-law exponent $b$), so it too survives regridding without a stored anchor. Both laws are shift-invariant: adding a constant to every $-\ln p$ leaves $k$ and the hard law's exponent $b$ unchanged.

`GridFn` is representation-agnostic: its explicit extrapolation API returns the literal $0.0$ for `Hard`. `DistributionFunction` owns the policies and maps this value to the policy-correct zero-probability value ($0$ under plain, $+\infty$ under negative-log). Every variant is therefore valid under both axes. The timetree passes store distributions on the negative-log axis (`Distribution<NegLog>`), where the ordinate is $-\ln p$ and normalization is a pure ordinate shift.

Builder methods `fn DistributionFunction.with_left_extrap()`, `fn DistributionFunction.with_right_extrap()`, and `fn DistributionFunction.with_extrap()` set the policy. Resampling preserves both policies. `fn DistributionFunction.negate_arg_inplace()` swaps left and right and transforms fitted laws because negating the argument reflects the domain.

`Distribution::Formula` evaluates its exact closure and stores no boundary policy. `fn Distribution.left_extrap()` and `fn Distribution.right_extrap()` therefore return `None` for Formula. Exact compact variants (`Empty`, `Point`, and `Range`) return `Some(Hard)`. The setters and `fn Distribution.fit_soft_tail()` change only Function distributions.

## Tail assignments in the timetree pipeline

The inference pipeline produces several distribution types. Each carries per-side tails that downstream operations use when computing support intersections. Tails are metadata set on the convolution result; they are consumed by the next operation in the pipeline (multiplication or division).

| Distribution            | Left tail (past)  | Right tail (future) | Where set                                                                                                 | Consumed by                                                                     |
| ----------------------- | ----------------- | ------------------- | --------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| Leaf time constraint    | `Error`           | `Error`             | default                                                                                                   | backward convolution input                                                      |
| Branch length           | `HardApproach`    | `Linear` (fitted)   | `fn compute_branch_length_distribution()` (hard approach near `t=0`, `SoftTailLaw::fit` near `t_max`)     | convolution input (negated for backward)                                        |
| Backward parent message | `Linear` (fitted) | `Hard`              | `fn Distribution::fit_soft_tail()` on the left, `with_right_extrap(Hard)` (backward_pass)                       | multiplication (combining children), division (forward cavity)                  |
| Internal node time dist | `Linear` (fitted) | `Hard`              | product result (left tail composed in closed form by `fn distribution_product()`, preserved by re-peaking) | forward convolution input, division dividend                                    |
| Forward message         | `Hard`            | `Linear` (fitted)   | `with_left_extrap(Hard)`, `fn Distribution::fit_soft_tail()` on the right (forward_pass)                        | multiplication (refining node dist)                                             |
| Refined node time dist  | `Hard` at the trim cut | `Hard` at the trim cut | multiplication result (composed from forward message and subtree dist, preserved by normalize), then cut at its mass quantiles ([timetree-forward-posterior-trimming.md](timetree-forward-posterior-trimming.md)) | `fn Distribution.likely_time()` extraction, forward cavity dividend |

### Why backward messages have `Linear` left / `Hard` right

The backward parent message is the convolution $\text{parent\_message} = \text{child\_time\_dist} \circledast (-\text{branch\_dist})$. It represents "when could the parent have existed, given this child?" The parent can be arbitrarily far in the past -- no child constrains how old its ancestor might be. The left side is therefore soft. It carries a fitted `Linear` tail: the message decays away from its peak, and a log-linear tail captures that decay with finite mass, whereas a flat tail would have infinite mass and no finite quantiles. The child's sampling date provides a hard upper bound on the parent's age (a parent cannot be more recent than its child), so the right tail is `Hard`.

### Why forward messages have `Hard` left / `Linear` right

The forward message is the convolution $\text{dist\_from\_parent} = \text{parent\_except\_subtree} \circledast \text{branch\_dist}$. It represents "when could this node have existed, given its parent and branch?" The parent's committed time provides a hard lower bound (branch lengths are non-negative), so the left tail is `Hard`. There is no upper bound from the parent side on how far in the future the node could be, so the right side is soft: a fitted `Linear` tail, for the same integrability reason as the backward left tail.

### Why the branch-length distribution has `HardApproach` left / `Linear` right

The branch-length distribution is the per-edge likelihood $L(t) \propto (\mu t)^{n} e^{-\mu t}$ over the branch duration, gridded on $[t_\text{min}, t_\text{max}]$. Its two edges are asymmetric facts. The left edge $t = 0$ is a hard boundary: a duration cannot be negative. For a branch with mutations the density vanishes there as a power law, so the left tail is a `HardApproach` law fitted near $t = 0$, which also carries the finite mode of a zero-mutation branch sitting on the boundary. The right edge $t_\text{max}$ is not a fact about the distribution, only where gridding stopped: the exponential $e^{-\mu t}$ keeps the density decaying past it (the Poisson indel term only adds a further linear neg-log slope), so the right side is a soft `Linear` tail fitted from the outermost points, with finite mass. Leaving it `Error` would declare a hard cutoff and drop any edge whose inferred duration exceeds $t_\text{max}$.

The FFT convolution of two functions reads only the grid ordinates and reconstructs the message tails from its own output, so it does not consume these tails ([kb/issues/M-distribution-convolution-drops-operand-tails.md](../issues/M-distribution-convolution-drops-operand-tails.md)); `negate()` swaps and sign-flips them for the backward message. The exact convolutions with a point or a range keep them (see [exact point and range convolution](#exact-point-and-range-convolution)). The tails enter the mass of the distribution (`fn Distribution.log_mass()`, the mass windows) and any evaluation beyond the grid: a duration beyond $t_\text{max}$ reads the extrapolated `Linear` tail instead of failing the `Error` default.

### Amplitude, re-peaking and normalization

A `DistributionFunction` stores a shape on its grid and a level in its amplitude, as `Point` and `Range` do: its values are `Y::multiply(ampl, stored)`. The amplitude defaults to the unit level ($0$ under negative-log), so a function built from values holds them as stored values.

- `fn Distribution.repeak()` moves the peak of the stored values into the amplitude and changes no value. It keeps the grid near the peak (no underflow) without losing the level. A function with no finite peak becomes `Empty`; a `NaN` value or amplitude is an error.
- `fn Distribution.normalize()` sets the peak to the unit level and drops the level.
- Under the negative-log axis both are an ordinate shift (`fn DistributionFunction.shift_y()`). A shift is closed form on both fitted laws -- the soft-tail slope and the hard approach law's shape are shift-invariant -- so the laws carry through unchanged. Every `DistributionFunction` regrid also preserves its per-side policy, fitted law and amplitude, and `fn DistributionFunction.scale_y()` rescales the fitted laws.
- Equality compares values: two functions with the same grid, tails and values are equal however the level is split between grid and amplitude.

The operations keep the amplitude:

| Operation | Result amplitude |
| --- | --- |
| convolution | product of the operands' amplitudes; the FFT runs on the stored values, the exact point and range convolutions add the amplitudes in neg-log |
| product and division of functions | product or quotient of the amplitudes; the stored values are combined on the common grid |
| `fn distribution_multiply_by_fn()` | the input's amplitude; the weights are added to the stored values |
| resampling, re-windowing to mass, negation, tail fitting | unchanged |

Combining stored values rather than values including the amplitude matters for reproducibility: adding a level of order one before re-peaking changes the last bits of the values, and the discrete grid choices downstream (point counts, mass windows) amplified that into node-date changes of up to $0.0025$ years.

The backward pass re-peaks each node's combined distribution, the edge convolution and the mass re-windowing re-peak their results, so every message and subtree distribution carries its likelihood as its level, coalescent factors included. `fn Distribution.log_mass()` integrates the shape with its tails and adds the amplitude; `fn log_marginal_likelihood()` evaluates it for the root's subtree distribution, the log marginal likelihood of the node times. Branch-length distributions stay peak-normalized, so the log marginal likelihood leaves out their peaks, which do not depend on the root or the clock rate. The forward pass normalizes posteriors, so posterior values are unchanged.

The log marginal likelihood is a quadrature of the backward messages, with these error sources:

- **Mass windows**: each backward message is resampled onto the window that drops $\epsilon = 5 \cdot 10^{-4}$ of its mass per side, and the soft tail law fitted at the window edge carries the cut mass. On the test fixtures, raising the message grid from 300 to 900 points and lowering $\epsilon$ to $10^{-9}$ changes no digit of the log marginal likelihood
- **Trapezoid quadrature**: the trapezoid rule of the grid samples has an error of second order in the spacing. Against exact piecewise integration of a three-tip tree with 2001-point branch distributions, the gap is $4.8 \cdot 10^{-6}$ to $8.2 \cdot 10^{-6}$, and $5.4 \cdot 10^{-7}$ to $9.2 \cdot 10^{-7}$ at 6001 points; a coalescent prior adds the second-order error of its kinks between grid points (a factor 4.09 per doubling). The tests of the marginal likelihood and of the clock rate step assert agreement with their closed-form oracles at $10^{-6}$ or tighter
- **Log-linear resampling**: products and FFT operands resample their operands log-linearly, which undercounts a Gaussian of standard deviation $\sigma$ by about $\Delta t^2 / (12\,\sigma^2)$ relative per resampling at the spacing $\Delta t$, an estimated $4 \cdot 10^{-5}$ on 300-point grids; this is not measured on real runs
- **Grid size of real runs**: with 300-point messages, the first-step log marginal likelihood on simulated data misses its converged value by 5 to 32 nats ([kb/issues/H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md](../issues/H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md)). Far below their peak the backward messages are extrapolated tail laws, which can be hundreds of nats too heavy ([kb/issues/H-timetree-backward-messages-too-heavy-far-below-their-peak.md](../issues/H-timetree-backward-messages-too-heavy-far-below-their-peak.md))

The backward pass combines child messages with a single product and re-peaking: the product composes the result tails (see [multiplication tails](#multiplication)) and re-peaking preserves them, so the folded result keeps its fitted `Linear` left tail and can still extend to reach a child with a disjoint finite grid. The child fold is `fn distribution_product()` (treetime-distribution), which co-locates the `Function` messages on one common grid and composes the summed left tail in closed form from the operands' own tails (`fn compose_multiplication_tail()`), so no tail is re-fit from the combined grid.

The forward pass applies its own tail policies _after_ re-peaking because the forward message needs `Hard` left and a fitted `Linear` right regardless of the tails the multiplication produced. It fits the right soft tail from the re-peaked grid; the shift left the slope unchanged, so fitting before or after re-peaking gives the same law.

## Grid intersection contract

v1 computes the result grid for both multiplication and division by resampling onto a new uniform grid spanning the exact analytical support intersection. This is an intentional divergence from v0, which collects knots from both operands inside the overlap and builds a non-uniform grid from them ([distribution.py#L82-L145](../../packages/legacy/treetime/treetime/distribution.py#L82-L145)).

Given an intersection $[x_{\min}, x_{\max}]$, the result grid has:

$$n = \max\left(\left\lceil\frac{x_{\max} - x_{\min}}{\Delta x}\right\rceil + 1,\; 2\right)$$

where $\Delta x = \min(\Delta x_a, \Delta x_b)$ for Function * Function, and $\Delta x = \Delta x_f$ for Range * Function or Formula * Function (using the Function operand's spacing).

The interval count is rounded **up** (the resolution floor): the realized spacing $(x_{\max} - x_{\min})/(n-1)$ is then never coarser than the finest operand's $\Delta x$, so multiplying a narrow distribution by a wide one cannot ratchet the spacing coarser than the narrow operand had. A ratio within a small relative tolerance of an integer is taken as that integer, so float error at an exact multiple does not add a spurious grid point.

The result grid includes both analytical endpoints via `Array1::linspace`, so a hard bound lands exactly on a grid node and a mode sitting on it stays representable. A disjoint intersection produces `Empty`. Endpoint-only contact produces a `Point` distribution evaluated at the shared boundary, matching v0's behavior of converting one surviving knot to a delta.

The same `Array1::linspace` grid and support trichotomy serve both the pairwise `fn multiply_function_function()` and the N-ary `fn multiply_functions()` (the shared co-location primitive that `fn distribution_product()` folds the child messages with), so a pairwise product equals the two-operand N-ary product.

Endpoint contact uses exact comparison of the analytical bounds, not a tolerance. A tolerance would enlarge the intersection and could turn two narrowly disjoint supports into a `Point`.

A count above the grid point limit stops the run with an error instead of an allocation ([distribution-grid-point-limit.md](distribution-grid-point-limit.md)).

### How tails extend the intersection

Both multiplication and division compute the intersection after extending each operand's evaluable domain based on its tail policy. The extension is per-side and independent: each side of each operand is resolved from its own `BoundaryBehavior`, and the result grid spans the intersection of the two extended domains. Without extension, operands with disjoint finite grids but overlapping via tails would collapse to `Empty`.

#### Multiplication

Each operand extends its own domain symmetrically. The extension limit is the other operand's grid bound on the same side:

| Operand tail | Behavior                                 | Rationale                                                                                                                                                |
| ------------ | ---------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Linear`     | Extend to the other operand's grid bound | The operand is evaluable at its exponential tail value beyond the grid; the product reflects the informative operand's shape scaled by the decaying tail |
| `Hard`       | Keep grid boundary                       | The value outside support is zero; $0 \cdot f(x) = 0$, so the product carries no information in the extension region                                     |
| `Error`      | Keep grid boundary                       | Out-of-support evaluation is undefined; strict finite-grid intersection                                                                                  |

Equivalently, per side the result takes the tightest (innermost) hard bound and the loosest (outermost) soft bound; a hard bound dominates a soft bound on the same side. Extending each soft side to the other operand's grid bound and intersecting selects this automatically.

The composed result tail per side is the strongest class of the two operands, `soft < hard < Error`, where the hard class covers both `Hard` and `HardApproach`. Two `Linear` tails compose in closed form: multiplication is addition in $-\ln p$, so their slopes add ($k = k_a + k_b$). Two `HardApproach` laws are not composed: no production path multiplies them (message hard sides are the nullary `Hard`, and only the branch-length factor produces `HardApproach`), and at distinct boundaries the product is not a single-exponent power law at the tighter bound, so that case is a loud internal error rather than a lossy composition. `Linear` always carries a fitted slope; the total enum has no unfitted `Linear`, so no composition step defers a refit.

When both operands have `Error` or `Hard` tails (the default), the result is identical to strict intersection. When operands carry soft `Linear` tails -- as backward parent messages do (see [tail assignments](#tail-assignments-in-the-timetree-pipeline)) -- disjoint finite grids overlap via the tails. This prevents the product from collapsing to `Empty` when temporal signals conflict, for example under `--keep-root` where rerooting cannot resolve the tension between subtrees.

Point * Function follows the same per-side rule without constructing a result grid. A point in a `Linear` tail evaluates against the tail's exponential value. A point beyond an `Error` or `Hard` boundary produces `Empty`.

##### Empty invariant

A Function-producing multiplication returns `Empty` only when the operands' hard domains are genuinely disjoint. The hard domain of an operand is its grid bounds on hard sides and unbounded ($\pm\infty$) on soft sides, since a soft side continues under its tail law rather than terminating. Because a soft side always bridges a gap to the other operand, an empty product can arise only from two hard bounds facing each other across a gap; any other empty result is a numerical or logic collapse that would silently poison every ancestor to the root (the original motivating defect). `fn multiplication_empty_result()` checks the hard domains and reports `fn make_internal_error!()` rather than returning `Empty` when they overlap.

##### Result tails

A `Function` result carries per-side tails composed from the two operands' tails on that side. Beyond a boundary the product is evaluated pointwise: if either operand is undefined there the product is undefined (`Error`); otherwise if either operand is zero the product is zero (a hard class); only when both operands continue softly is the product soft. This is the maximum over the precedence `soft` < `hard` < `Error`, where soft covers `Linear` and hard covers `Hard`/`HardApproach` (the more restrictive tail wins):

| A tail  | B tail  | Result tail                                            |
| ------- | ------- | ------------------------------------------------------ |
| soft    | soft    | soft (`Linear`, slopes added)                          |
| soft    | hard    | hard                                                   |
| soft    | `Error` | `Error`                                                |
| hard    | hard    | hard                                                   |
| hard    | `Error` | `Error`                                                |
| `Error` | `Error` | `Error`                                                |

A `Range` operand has exact hard boundaries. A `Formula` operand has no stored policy; arithmetic uses its finite sampling range as an `Error`-bounded domain when it must produce a Function result. Range * Function and Formula * Function results therefore do not inherit a soft tail from the compact or Formula operand. Only Function * Function can produce soft result tails, and only when both operands opt in. Combined with tail-preserving normalization, this is what lets the backward pass accumulate child messages without re-applying tails.

#### Division

Division shares multiplication's support rule. The forward pass divides a parent's time distribution (dividend $f$) by a child's backward message (divisor $g$) to compute the cavity $f/g$. In neg-log the quotient is a subtraction of ordinates ($-\ln(f/g) = -\ln f - (-\ln g)$).

The dividend is the product of _all_ the parent's factors, and the divisor is _one_ of those factors: $f = g \cdot \text{rest}$, so the cavity $f/g = \text{rest}$ is itself a product of factors. Because the divisor is a factor of the dividend, its decay is already baked into the dividend and cancels in the quotient. On a soft side the neg-log slopes add under multiplication (see [multiplication tails](#multiplication)), so $k_f = k_g + k_\text{rest}$ and the quotient tail slope is $k_f - k_g = k_\text{rest} \ge 0$: the quotient decays, it does not spike. Dividing by $g$ _in isolation_ would inflate the result where $g$ decays, but the cavity never does that -- $g$ appears in $f$.

The result grid is therefore the same per-side rule as multiplication (intersect hard, union soft):

- any hard operand on a side (dividend **or** divisor) bounds the quotient at the innermost hard edge;
- otherwise both operands are soft, and the quotient extends to the outermost soft edge, sampling each operand through its own tail law and subtracting.

A soft divisor edge no longer truncates the quotient: the divisor's decaying tail is evaluated as bulk out to the union edge, exactly as the Python prototype `test_scripts/density_algebra.py` `combine()` does with `signs=[1, -1]`. A divisor `Hard`/`HardApproach`/`Error` edge strictly inside the dividend still bounds the quotient there -- dividing by zero (`+\infty` under neg-log, giving $f - (+\infty) = -\infty$) or an undefined value is a genuine spike -- and the quotient is `Error` beyond it.

##### Result tails

Per side the quotient tail follows from which operand binds the edge:

| Dividend side | Divisor side | Binding edge   | Result tail                                                                                                  |
| ------------- | ------------ | -------------- | ------------------------------------------------------------------------------------------------------------ |
| soft          | soft         | outermost soft | refit `Linear` from the combined grid                                                                        |
| hard          | soft         | dividend hard  | the dividend's own hard tail (`Hard`/`HardApproach`), or `Error` if the dividend side is `Error`             |
| soft          | hard         | divisor hard   | `Error` (undefined or a spike past the divisor edge)                                                         |
| hard          | hard         | innermost hard | the dividend tail if the dividend binds (ties included, a $0/0$ cavity edge the dividend owns), else `Error` |

A soft result side is **refit** from the combined grid via `fn SoftTailLaw::fit()`, not composed in closed form. Exact slope subtraction $k_f - k_g$ holds only while every operand's grid reaches the same edge; the cavity's operands are independently regridded (each fits its own tail from FFT output), so the outermost quotient cells mix one operand's bulk with another's tail and the true local slope there is not $k_f - k_g$. The refit reads the quotient's own outermost points, and `fn SoftTailLaw::fit()` clamps a wrong-sign slope to flat, which is the guard against a tail that would manufacture probability outward. Multiplication composes its soft tails in closed form (slopes add, always decaying at least as fast) because addition cannot invert a decay; division refits because subtraction can.

The soft sides union and extend under their tail laws, so an independently re-gridded dividend (for example, a mass-sized parent posterior whose grid extends past a hard-bounded child message) keeps the cavity non-empty even when the two finite grids barely overlap.

A subtlety this rule does not correct: when the divided-out child held the _tightest_ hard bound on a side, the parent product never stored the next factor's bound, so the quotient can only keep the dividend's (former child's) hard bound there -- a $0/0$ cancellation at the edge. The cavity is then slightly too tight on that side. Restoring the exact bound would require rebuilding the cavity from the retained factors instead of dividing the aggregate, tracked in [kb/issues/N-timetree-cavity-hard-bound-too-tight.md](../issues/N-timetree-cavity-hard-bound-too-tight.md).

Coalescent contributions do not introduce another grid. The backward pass multiplies child messages first and evaluates the coalescent contribution pointwise on the resulting grid in negative-log space (`fn distribution_multiply_by_fn()`).

### Tails of a weighted distribution

`fn distribution_multiply_by_fn()` [packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs](../../packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs) composes each soft tail law with the weight: the new neg-log slope is the incoming slope plus the slope of the weight just beyond the grid edge. The owner of the weight states that slope from its own extrapolation rule (the `outward_slope` argument); for the coalescent factor it is `fn CoalescentPrior::node_contribution_outward_slope()` [packages/treetime/src/coalescent/mean_field.rs#L125](../../packages/treetime/src/coalescent/mean_field.rs#L125). A composed slope that does not decay is an internal error that names the side and both slopes; a probe on `data/flu/h3n2/200` and `data/sc2/4500` with `--coalescent-opt` and on the coalescent runs of the simulated validation subset found none in 96 290 compositions (largest composed slope $-0.41$). Hard sides keep their law.

The slope is not taken from the grid points. The coalescent cost steps wherever the number of lineages changes, and going into the past it can step down, where an older sample adds a lineage. A refit from the outermost grid points that straddles such a step fits a slope rising outward, which the clamp of `fn SoftTailLaw::fit()` turns flat: the node distribution then has infinite mass. This happened at 3-7 recent nodes per 400-tip simulation with a coalescent prior. A secant of sampled weights fails the same way: a bounded step inside the sampled span gives a slope that does not decay, while the product has finite mass.

## Mass of a distribution

`struct MassProfile` [packages/treetime-distribution/src/distribution_ops/mass_profile.rs](../../packages/treetime-distribution/src/distribution_ops/mass_profile.rs) integrates a `DistributionFunction` with its tails: regular cells by the trapezoid rule of the grid samples, a cell with one zero-density end point by the power law that evaluation uses in that cell (`fn GridFn::hard_cell()`), and the `Linear` and `HardApproach` tail laws in closed form. `fn Distribution.log_mass()`, the mass windows, the forward-posterior trim, the root median, the survival functions and the 90% intervals all read it ([timetree-posterior-statistics-count-tail-mass.md](timetree-posterior-statistics-count-tail-mass.md)).

A mass window (`fn mass_bounded_domain()`) moves only a `Linear` side inward, to the quantile on the left or the inverse survival on the right; a `Hard` or `HardApproach` side keeps its grid edge whatever its density, and an `Error` side is an error. The forward-posterior trim is the only operation that moves a hard side inward.

## Exact point and range convolution

A date range $[a, b]$ convolved with a function $f$ is computed exactly from the mass model of $f$ (`fn convolution_range_function()` [packages/treetime-distribution/src/distribution_ops/convolve.rs#L79-L137](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L79-L137)):

$$g(t) = A_r\,A_f \int_{t - b}^{t - a} f(s)\,ds$$

where:

- $g$ -- density of the result
- $A_r$, $A_f$ -- amplitudes of the range and of the function

The result spans the full support, from the left end of $f$ plus $a$ to the right end of $f$ plus $b$, where an end is the grid edge for `Linear`, `Hard` and `Error` and the hard bound for `HardApproach`. Its grid has at most the spacing of $f$; when hard cells occur (a physical end, or a zero-density gap), the spacing is at most $1/8$ of the smaller of $b - a$ and the narrowest support component of $f$, so the rise of the density from each end is resolved. The physical end points stay on the grid with zero density, so the end cells are hard cells whose power-law mass carries the density that rises from each end. The result keeps the tail laws of $f$ on its soft sides and carries the summed amplitude, which the log marginal likelihood integrates.

A point convolution shifts the function, its tail laws and the bound of a `HardApproach` law, and adds the amplitudes. A range convolved with a range is an error. Point and range results skip the mass window of `fn convolve_across_edge()`, which would cut their exact soft law; every result's non-soft side is `Hard` at its grid edge.

A shift of the function to the midpoint of the range on the function's own grid, the earlier v1 rule, cuts half the range width from each side, so a parent is capped at the midpoint of a tip's date range. Against that rule, in the smoke matrix, three year-only tips of `rsv/a/100` move to the end of their range by 1.0 to 1.23 of their 90% width, and every other changed date stays within its 90% width. On `data/dengue/500` (325 of 457 tips with year-only dates) the runs take 6.7% more grid points and the largest grid grows from 670 to 2365 points, with dates moving by up to 0.137 years; CPU time falls on `data/dengue/500` (147.6 s to 104.2 s) and `data/mpox/clade-ii/500` (9.80 s to 6.49 s). Under the midpoint rule, a leaf range message with the left slope $-4.6 \cdot 10^{-12}$ put the $\epsilon$ quantile at $-1.66 \cdot 10^{12}$ years on `data/mpox/clade-ii/1000` at a clock rate of $10^{-3}$; with the exact range convolution that run completes.

## Divergence from v0

Two specific divergences from v0's distribution handling:

1. **Convolution tail reconstruction.** v0 rebuilds convolution tails inside `fn NodeInterpolator.convolve_fft()` ([node_interpolator.py#L231-L256](../../packages/legacy/treetime/treetime/node_interpolator.py#L231-L256)) as slope-based extrapolation (linear in neg-log, exponential in probability) from the outermost trusted points. v1 reconstructs the same tails in `fn convolution_function_function_fine()`: the FFT runs in plain probability space and is trusted only above a peak-relative floor (`1e-13`, matching v0), and the sub-floor tail is rebuilt by a two-point secant slope in negative-log space, extended only where it decays away from support. Each operand is converted to peak-normalized plain probability for the FFT (`fn to_peak_normalized_plain()`) and the result back to negative log. The inference passes then fit a per-message soft tail from the reconstructed grid (`fn Distribution::fit_soft_tail()`) and attach it as the message's out-of-support policy on the soft side, keeping the hard side `Hard`. That policy governs how the _following_ multiplication or division extends the support intersection.

2. **Default out-of-support handling.** v0 returns effectively zero probability outside support as a soft value (`fill_value=1e10` in neg-log, i.e. $\exp(-10^{10})$). v1 returns `Error` by default, requiring an explicit tail opt-in.

The v1 pass soft tails serve two ends. They preserve node-time monotonicity -- without a soft left tail on the backward message a parent's inferred time can be truncated too recent, producing negative branch lengths, and the soft tail lets the combined distribution place the parent appropriately older. The fitted `Linear` tail also has finite mass, so the quantile and interval integrals stay finite; because log-linear extrapolation slightly over-estimates the far tail it can never manufacture a spurious zero, so the monotonicity guarantee is retained. A separate forward-pass projection currently clamps committed internal-node point estimates; the statistical contract for that projection remains open in [kb/issues/M-timetree-marginal-node-times-can-violate-topology.md](../issues/M-timetree-marginal-node-times-can-violate-topology.md).

The hard side of each message uses the nullary `Hard` variant, which places zero probability past the grid edge. A mode sitting exactly on a hard boundary (as for a zero-mutation branch) is carried by placing the boundary as an exact grid endpoint, not by a fitted law: the branch-length factor stores a finite `t = 0` boundary as an exact endpoint with a `Hard` left side, so its mode lands on a real grid point.

## Accepted limitations

This is a bounded resampling rule, not a derived interpolation-error criterion. Rounding the interval count up keeps the realized spacing at or below the finest operand's $\Delta x$. Sequential (chained pairwise) multiplication still resamples at each binary operation and is not guaranteed to be associative. The N-ary `fn multiply_functions()` fold, by contrast, co-locates all operands on one grid in a canonical operand order, so its result is bit-identical across any permutation of the factors.

The FFT convolution runs in plain probability on peak-normalized operands and converts back to negative log, so its result carries ULP-level round-trip noise from the `-ln`/`exp` bracketing, far below any scientific tolerance.

Scientific workflows requiring posterior peak, normalization, or integrated-probability error bounds need a separate accuracy contract.

## Code references

### Tail policy

- [packages/treetime-grid/src/boundary_behavior.rs](../../packages/treetime-grid/src/boundary_behavior.rs): `enum BoundaryBehavior` (`Error`/`Hard`/`HardApproach`/`Linear`, every variant a complete value), `fn BoundaryBehavior.is_soft()`/`is_hard()`, `const DEFAULT_TAIL_FIT_POINTS`
- [packages/treetime-grid/src/soft_tail_law.rs](../../packages/treetime-grid/src/soft_tail_law.rs): `struct SoftTailLaw` (single neg-log slope), `fn SoftTailLaw::fit()`, `fn SoftTailLaw::mass()`, `fn SoftTailLaw::compose_multiply()`
- [packages/treetime-grid/src/grid_fn.rs](../../packages/treetime-grid/src/grid_fn.rs): policy-free grid interpolation; `fn GridFn.interp_with_extrap()` and `fn GridFn.resample_with_extrap()` accept explicit policies for one operation
- [packages/treetime-distribution/src/distribution_core/function.rs](../../packages/treetime-distribution/src/distribution_core/function.rs): `left_extrap`/`right_extrap` fields and getters; `ampl` with `fn DistributionFunction.amplitude()`, `fn DistributionFunction.with_amplitude()`, `fn DistributionFunction.values()`, `fn DistributionFunction.repeak()`, value-based `PartialEq`; `fn DistributionFunction.interp()` (values) and `fn DistributionFunction.interp_stored()`; policy builders; tail fitting; policy- and amplitude-preserving resampling, ordinate transforms, and argument negation
- [packages/treetime-distribution/src/distribution_core/formula.rs](../../packages/treetime-distribution/src/distribution_core/formula.rs): exact Formula closure and finite sampling range, with no stored boundary policy
- [packages/treetime-distribution/src/distribution_core/distribution.rs](../../packages/treetime-distribution/src/distribution_core/distribution.rs): optional boundary getters; Function-only policy builders and tail fitting; hard compact variants; `fn Distribution.normalize()`, `fn Distribution.repeak()`, `fn Distribution.log_mass()`
- [packages/treetime-distribution/src/distribution_ops/mass_profile.rs](../../packages/treetime-distribution/src/distribution_ops/mass_profile.rs): `struct MassProfile`, `struct MassFractions`, `enum TailMass`
- [packages/treetime-distribution/src/distribution_ops/mass_domain.rs](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs): `fn log_mass()` (shape and tails plus amplitude), `fn mass_bounded_domain()`, `fn trim_with_profile()`
- [packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs](../../packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs): `fn distribution_multiply_by_fn()` composes the soft tail laws with the outward slope of the weight and keeps the amplitude
- [packages/treetime-distribution/src/distribution_ops/convolve.rs](../../packages/treetime-distribution/src/distribution_ops/convolve.rs): `fn convolution_range_function()`, `fn convolution_point_function()`, `fn convolution_function_function_fine()`

### Tail-aware arithmetic

- [packages/treetime-distribution/src/distribution_ops/multiply.rs](../../packages/treetime-distribution/src/distribution_ops/multiply.rs): `fn multiplication_support_intersection()`, `fn multiplication_empty_result()` (empty invariant), `fn hard_domains_disjoint()`, `fn compose_multiplication_tail()`, `fn with_composed_tails()`, `fn multiply_point_function()`, `fn multiply_function_function()`, `fn multiply_range_function()`
- [packages/treetime-distribution/src/distribution_ops/divide.rs](../../packages/treetime-distribution/src/distribution_ops/divide.rs): `fn divide_function_by_function()`, `fn divide_range_by_function()`, `fn divide_point_by_function()`, `fn apply_division_tail()`, `fn division_side_tail()` (reuses `fn multiplication_support_intersection()` for the grid, refits soft result tails)

### Inference pass tail application

- [packages/treetime-distribution/src/distribution_core/distribution.rs#L173](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L173): `fn Distribution::fit_soft_tail()` -- fits the soft-side `Linear` tail of a message, errors on a degenerate grid, inert on non-Function messages. Shared by both passes.
- [packages/treetime/src/timetree/inference/backward_pass.rs](../../packages/treetime/src/timetree/inference/backward_pass.rs): backward message soft-left / hard-right assignment; the child fold delegates to `fn distribution_product()` (treetime-distribution), which composes the summed tails in closed form; node distributions are re-peaked.
- [packages/treetime/src/timetree/inference/result.rs](../../packages/treetime/src/timetree/inference/result.rs): `fn log_marginal_likelihood()`
- [packages/treetime-distribution/src/distribution_ops/multiply.rs](../../packages/treetime-distribution/src/distribution_ops/multiply.rs): `fn distribution_product()` folds the child `Function` messages on one common grid; `fn compose_multiplication_tail()` composes each side's result tail in closed form (soft slopes add, a hard bound dominates).
- [packages/treetime/src/timetree/inference/forward_pass.rs](../../packages/treetime/src/timetree/inference/forward_pass.rs): forward message hard-left / soft-right assignment; posteriors are normalized, including the subtree distributions that the root and exactly dated nodes take as their posterior.

### Tests

- [`packages/treetime-distribution/src/distribution_ops/__tests__/test_multiply.rs`](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_multiply.rs): tail combinations (overlapping, disjoint, mixed `Linear`/`Hard`/`Error`, commutativity)
- [`packages/treetime-distribution/src/distribution_ops/__tests__/test_divide.rs`](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_divide.rs): exact intersection under default boundaries, extension under explicit divisor tails
- [`packages/treetime-distribution/src/distribution_core/__tests__/test_boundary_behavior.rs`](../../packages/treetime-distribution/src/distribution_core/__tests__/test_boundary_behavior.rs): boundary-variant evaluation and composition
- [`packages/treetime-distribution/src/distribution_ops/__tests__/test_amplitude.rs`](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_amplitude.rs): re-peaking keeps values, normalization drops the level, convolution adds log masses, multiplication adds values, `fn Distribution.log_mass()`
- [`packages/treetime-distribution/src/distribution_ops/__tests__/test_multiply_by_fn.rs`](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_multiply_by_fn.rs): the slope of the weight added beyond the soft edge, a step one cell inside the edge keeps a decaying tail, `fn Distribution.log_mass()` counts the composed tail
- [`packages/treetime-distribution/src/distribution_ops/__tests__/test_convolve.rs`](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_convolve.rs): range convolution against the box convolution of the trapezoid model and at second order against box times exponential and box times Gaussian, full support, both amplitudes, end cells
- [`packages/treetime-distribution/src/distribution_ops/__tests__/test_mass_profile.rs`](../../packages/treetime-distribution/src/distribution_ops/__tests__/test_mass_profile.rs) and the other `test_mass_profile*` files: cells, hard cells and tails against closed forms, convergence, properties
- [`packages/treetime/src/timetree/inference/__tests__/test_backward_pass.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_backward_pass.rs): `test_backward_pass_root_carries_the_marginal_likelihood` against the closed form for two dated leaves with Gaussian branch likelihoods
- [`packages/treetime-distribution/src/distribution_core/__tests__/test_function.rs`](../../packages/treetime-distribution/src/distribution_core/__tests__/test_function.rs): `fn Distribution::fit_soft_tail()` -- slope recovery, finite tail mass, inert on non-Function, error on degenerate grid

## Related knowledge base items

### Specifications

- [kb/algo/timetree-extrapolation-and-time-clamping.md](../algo/timetree-extrapolation-and-time-clamping.md): original spec defining extrapolation policies, exact intersection boundaries, and the monotonicity clamp. Source material for this decision.

### Proposals

- [kb/proposals/exponential-branch-length-tails.md](../proposals/exponential-branch-length-tails.md): a `BoundaryBehavior::Exponential` variant for branch length distributions, continuing the log-linear slope at the grid boundary.

### Issues

- [kb/issues/M-timetree-marginal-node-times-can-violate-topology.md](../issues/M-timetree-marginal-node-times-can-violate-topology.md): independent marginal modes can invert parent-child ordering. The forward-pass monotonicity clamp addresses this but its statistical contract is unresolved.
- [kb/issues/N-timetree-gm-runner-dataset-cases-disabled.md](../issues/N-timetree-gm-runner-dataset-cases-disabled.md): golden-master runner cases of several datasets stay disabled; tail-aware multiplication removed their original cause (empty distributions from disjoint multiplication).
- [kb/issues/N-distribution-mixed-nan-policy-undecided.md](../issues/N-distribution-mixed-nan-policy-undecided.md): `fn Distribution.normalize()` returns `Empty` for non-positive or non-finite max values, silently discarding error information.
