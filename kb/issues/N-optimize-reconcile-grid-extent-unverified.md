# `reconcile_zero_boundary` grid extent argument is unverified

> [!IMPORTANT]
> **Investigation required.** No test proves that the grid covers the true mode when the input branch length is uninformative. Build a fixture with the input branch length $100\times$ away from the true mode in each direction, as in step 1 of "Proposed action". The Dinh and Matsen K80 surface already fails the upper direction (see Evidence). Not measured: whether real data produces multimodal edges on `mpox/clade-ii/100`, `dengue/100`, datasets above 200 tips, and fixed `k80`/`hky85` models.

## Problem

`fn reconcile_zero_boundary()` [packages/treetime/src/optimize/zero_boundary.rs#L42](../../packages/treetime/src/optimize/zero_boundary.rs#L42) takes a `branch_length_for_extent` parameter separate from `candidate`. The separation exists because `candidate` may be clamped to $0$ or to a tiny floor like $10^{-12}$ after the inner solver runs, and the grid extent must remain tied to a meaningful scale. `fn optimize_edge()` [packages/treetime/src/optimize/dispatch.rs#L303-L304](../../packages/treetime/src/optimize/dispatch.rs#L303-L304) passes the solver start as the extent: the input branch length, or `one_mutation` or $k / \text{rate}$ when the input is zero and the edge needs a positive start.

This choice is correct in the common case: the input branch length is the optimizer's best prior estimate of the true scale. The corner case where the input is uninformative or wildly misleading is not covered:

- Input branch length is zero (e.g. uncalibrated tree): the grid spans $[0.1 \cdot \text{one\_mutation},\ \max(\text{one\_mutation}, 0.5)]$, which may not contain the true mode if the true scale is outside that interval.
- Input branch length is much smaller than the true mode (e.g. bad initial guess): `1.5 * bl + one_mutation < true_mode`, and the cap at `0.5` saves the day only if `true_mode < 0.5`.
- Input branch length is much larger than the true mode: the grid lower end stays at $0.1 \cdot \text{one\_mutation}$, so the 100 points spread over more decades and the resolution per decade drops; the true mode is missed only when it lies below the lower end.

## Evidence

Throwaway instrumentation of the per-edge optimizer on `rust` at `77dd3a84`, described in [M-optimize-brent-setup-differs-from-v0.md](M-optimize-brent-setup-differs-from-v0.md#evidence):

- **The extent misses the global maximum on the Dinh and Matsen K80 test surface.** The surface used in `test_dispatch_zero_boundary.rs` has a local maximum at $t \approx 0.167$ and its global maximum at $t \approx 6.37$. For every start $t_0 \leq 2$, the grid upper end $\max(1.5 t_0 + \text{one\_mutation}, 0.5)$ is below $6.37$, so both the grid search and the Brent bracket return the local maximum or the grid's own upper end ($t = 3.01$ from $t_0 = 2$). The existing tests assert the local maximum near $0.2$ as the expected result
- **The grid lower end misses maxima below $0.1 \cdot \text{one\_mutation}$.** On 4 Newton optimizations of `ebola/100` (start $t = 0$), the maximum lies at $t \approx 2.4\cdot10^{-6}$, below the grid lower end $5.3\cdot10^{-6}$, and the grid search returns zero with a log-likelihood shortfall of $1.5\cdot10^{-3}$
- **On the measured example datasets no edge had more than one local maximum** (`flu/h3n2/{20,200}`, `ebola/{20,100}`, `zika/86`, `rsv/a/100`; dense and sparse; `--model infer|jc69`; both methods). The extent question therefore matters for non-unimodal surfaces and very short branches, not for the common case on these datasets

> [!WARNING]
> **Needs review.** The 4 `ebola/100` cases came from the grid fallback of the Newton method, which `fn newton_sqrt_maximize()` [packages/treetime/src/optimize/method_newton.rs](../../packages/treetime/src/optimize/method_newton.rs) removed. The grid lower end still applies to the grid of `fn reconcile_zero_boundary()`

## Impact

Negligible when the initial guess is calibrated (the standard workflow). On datasets with uncalibrated input trees or bad initial guesses combined with non-unimodal models that force `reconcile_zero_boundary` to grid-search, the reconciled branch length lands on whichever grid point beats zero, and there is no test that proves the true mode is on the grid. The behavior is unverified.

## Proposed action

1. Test with a fixture where the input branch length is $100\times$ away from the true mode in each direction (above and below), verify the grid still covers the true mode.
2. If the test fails, consider passing a more informative extent: e.g. the maximum of the input branch length and the standard deviation of the Hamming distance, $\sqrt{p (1 - p) / L}$ for the observed difference fraction $p$ over $L$ sites (at most $1 / (2 \sqrt{L})$), or an invariant derived from the current contribution coefficients.
3. Document the admissibility contract: "the grid covers $[0.1 \cdot \text{one\_mutation},\ \max(1.5 \cdot \text{extent} + \text{one\_mutation}, 0.5)]$; the true optimum must lie inside it".

## Cross-references

- [packages/treetime/src/optimize/zero_boundary.rs](../../packages/treetime/src/optimize/zero_boundary.rs) (`reconcile_zero_boundary`, `grid_search_branch_lengths`)
- [N-optimize-zero-boundary-handled-in-several-places.md](N-optimize-zero-boundary-handled-in-several-places.md): option 2 there needs the extent decided here
- [N-optimize-grid-search-upper-bound-capped.md](N-optimize-grid-search-upper-bound-capped.md)

## Readiness

The proposed test is evidence gathering, and the alternative extent formulas are speculative. Start implementation only after a failing fixture identifies the required admissible interval.
