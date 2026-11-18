# Grid search resolution (100 points) is unverified

## Problem

`fn grid_search_branch_length()` [packages/treetime/src/optimize/zero_boundary.rs#L79-L108](../../packages/treetime/src/optimize/zero_boundary.rs#L79-L108) enumerates `GRID_SEARCH_POINTS = 100` [zero_boundary.rs#L9](../../packages/treetime/src/optimize/zero_boundary.rs#L9) log-spaced positive candidates between `0.1 * one_mutation` and `fn search_upper_bound()` $= \max(1.5 \cdot \text{start} + \text{one\_mutation}, 0.5)$ (`fn grid_search_branch_lengths()` [zero_boundary.rs#L110-L122](../../packages/treetime/src/optimize/zero_boundary.rs#L110-L122)). It returns the best grid point, or zero when zero is better, without local refinement. The only caller is `fn reconcile_zero_boundary()` [zero_boundary.rs#L57-L77](../../packages/treetime/src/optimize/zero_boundary.rs#L57-L77), for both the Brent and the Newton method. The 100-point choice is neither derived from a convergence analysis nor tested for sufficiency. A multi-modal surface with a narrow mode could fall between adjacent grid points and get missed, and the reconciliation path in `reconcile_zero_boundary` could then prefer a suboptimal grid point over a true positive mode.

The log spacing spans 3-4 decades for typical branch lengths (from $\sim 10^{-4}$ to $\sim 0.5$), so 100 points give roughly 25-33 points per decade. A mode narrower than one grid step would be invisible.

## Evidence

- **Spacing.** Consecutive grid points differ by the factor $(\text{upper} / \text{lower})^{1/99}$: 1.10 to 1.12 on the example datasets `flu/h3n2/{20,200}`, `ebola/{20,100}`, `zika/86`, `rsv/a/100`, where the grid spans $[0.1 \cdot \text{one\_mutation}, 0.5]$. The best grid point can therefore lie about 5% away from a maximum even for a broad mode. A short local refinement (Brent between the neighbors of the best grid point) removes this error without more grid points
- **Grid results on the reconcile path.** Throwaway instrumentation on `rust` at `77dd3a84` (method in [M-optimize-brent-setup-differs-from-v0.md](M-optimize-brent-setup-differs-from-v0.md#evidence)) recorded 752 Brent optimizations where the reconcile step changed the solver output, all to zero (same datasets, dense and sparse, `--model infer|jc69`). On these runs the grid never returned a positive point, so its resolution did not affect a positive result

> [!IMPORTANT]
> **Investigation required.** The Brent figures above were not re-measured on the current code. Not measured: how often the reconcile grid search returns a positive point, how far that point lies below the reference maximum, and the same for the current safeguarded Newton method, on the datasets above plus `mpox/clade-ii/100`, `dengue/100`, and fixed `k80`/`hky85` models.

## Impact

Negligible for surfaces with modes wider than one grid step. Unknown for surfaces with sharp peaks. The Dinh and Matsen 2017 K80 counterexample has a broad local max near $t \approx 0.2$ which is well-covered, but non-unimodal GTR models with extreme parameter values could produce sharper modes.

## Proposed action

1. Convergence test: for a fixture surface with a known mode (e.g. a Gaussian bump of width $w$), vary `GRID_SEARCH_POINTS` across $\{20, 50, 100, 200, 500, 1000\}$ and measure the grid's mode-detection accuracy. Establish the minimum grid density that catches a mode of width $w$.
2. Document the result as "grid search catches modes of width $\geq X$ subs/site per decade".
3. If 100 points is insufficient for realistic sharp modes, either raise the constant or switch to a two-pass grid (coarse scan to locate a concave region, then refinement around the best candidate).

## Cross-references

- [packages/treetime/src/optimize/zero_boundary.rs](../../packages/treetime/src/optimize/zero_boundary.rs) (`grid_search_branch_length`, `grid_search_branch_lengths`, `GRID_SEARCH_POINTS`)

## Proposed implementation plan

The locations above are current: `GRID_SEARCH_POINTS = 100` at `zero_boundary.rs#L9`, no refinement after `grid_search_branch_length()`.

1. Measure first (the Investigation required block): add temporary instrumentation (stderr) to `reconcile_zero_boundary()` and count, per dataset and method, how often the grid returns a positive point and its relative distance to a dense reference maximum (a 10 000-point log grid plus Brent on the same objective). Remove the instrumentation afterwards.
2. If the grid never returns a positive point on the listed datasets, record the result in this issue and close it as negligible. Otherwise continue.
3. Add a local refinement in `grid_search_branch_length()`: Brent (`argmin`, as the main Brent path uses) on the bracket of the two neighbors of the best grid point, and keep the better of the refined point and zero.
4. Tests in `packages/treetime/src/optimize/__tests__/test_grid_search/`: (a) a parameterized `rstest` on synthetic objectives with a known maximum (a Gaussian bump in log space of widths 0.5, 0.1, and 0.02 decades, maximum known analytically), asserting the returned point within `1e-8` relative of the analytic maximum with refinement; (b) the narrowest width that the 100-point grid catches, as a documented boundary case. Oracle: the analytic maximum of the synthetic objective.
5. Validate with `just test-rs`, and `just smoke` only if outputs change on the datasets (per `cli-verification.instructions.md`).
