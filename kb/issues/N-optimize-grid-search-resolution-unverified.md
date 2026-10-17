# Grid search resolution (100 points) is unverified

## Problem

`fn grid_search_branch_length()` [packages/treetime/src/optimize/zero_boundary.rs#L67](../../packages/treetime/src/optimize/zero_boundary.rs#L67) enumerates `GRID_SEARCH_POINTS = 100` log-spaced positive candidates between `0.1 * one_mutation` and `max(1.5 * bl + one_mutation, 0.5)`. The 100-point choice is neither derived from a convergence analysis nor tested for sufficiency. A multi-modal surface with a narrow mode could fall between adjacent grid points and get missed, and the reconciliation path in `reconcile_zero_boundary` could then prefer a suboptimal grid point over a true positive mode.

The log spacing spans 3-4 decades for typical branch lengths (from $\sim 10^{-4}$ to $\sim 0.5$), so 100 points give roughly 25-33 points per decade. A mode narrower than one grid step would be invisible.

## Evidence

The grid result is used without refinement. Throwaway instrumentation on `rust` at `77dd3a84` (method in [M-optimize-brent-setup-differs-from-v0.md](M-optimize-brent-setup-differs-from-v0.md#evidence)) measured the Newton path that returns the grid result directly when the start is not concave:

- 17 of 9,544 Newton optimizations took this path (`flu/h3n2/{20,200}`, `ebola/{20,100}`, `zika/86`, `rsv/a/100`; dense and sparse; `--model infer|jc69`)
- 15 of them ended more than $10^{-6}$ log-likelihood units below the reference maximum, 3 more than $10^{-3}$. Example: `flu/h3n2/200`, grid point $7.28\cdot10^{-4}$, reference $7.11\cdot10^{-4}$, shortfall $2.7\cdot10^{-4}$
- Consecutive grid points differ by a factor of 1.10 to 1.12 on these datasets (the grid spans $[0.1 \cdot \text{one\_mutation}, 0.5]$ in 100 points), so the best grid point can be up to about 5% away from a maximum even for a broad mode A short local refinement (Brent between the neighbors of the best grid point) removes this error without more grid points

> [!WARNING]
> **Needs review.** The Newton path measured here no longer exists: `fn newton_sqrt_maximize()` [packages/treetime/src/optimize/method_newton.rs](../../packages/treetime/src/optimize/method_newton.rs) replaced the grid fallback with a bracketed search. The grid now runs only inside `fn reconcile_zero_boundary()`, for both methods

## Impact

Negligible for surfaces with modes wider than one grid step. Unknown for surfaces with sharp peaks. The Dinh and Matsen 2017 K80 counterexample has a broad local max near $t \approx 0.2$ which is well-covered, but non-unimodal GTR models with extreme parameter values could produce sharper modes.

## Proposed action

1. Convergence test: for a fixture surface with a known mode (e.g. a Gaussian bump of width $w$), vary `GRID_SEARCH_POINTS` across $\{20, 50, 100, 200, 500, 1000\}$ and measure the grid's mode-detection accuracy. Establish the minimum grid density that catches a mode of width $w$.
2. Document the result as "grid search catches modes of width $\geq X$ subs/site per decade".
3. If 100 points is insufficient for realistic sharp modes, either raise the constant or switch to a two-pass grid (coarse scan to locate a concave region, then refinement around the best candidate).

## Cross-references

- [packages/treetime/src/optimize/zero_boundary.rs](../../packages/treetime/src/optimize/zero_boundary.rs) (`grid_search_branch_length`, `grid_search_branch_lengths`, `GRID_SEARCH_POINTS`)
