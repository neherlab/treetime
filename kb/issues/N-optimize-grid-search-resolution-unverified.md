# Grid search resolution (100 points) is unverified

## Problem

`grid_search_inner` at [packages/treetime/src/commands/optimize/optimize_unified.rs](../../packages/treetime/src/commands/optimize/optimize_unified.rs) enumerates `GRID_SEARCH_POINTS = 100` log-spaced positive candidates between `0.1 * one_mutation` and `max(1.5 * bl + one_mutation, 0.5)`. The 100-point choice is neither derived from a convergence analysis nor tested for sufficiency. A multi-modal surface with a narrow mode could fall between adjacent grid points and get missed, and the reconciliation path in `reconcile_zero_boundary` could then prefer a suboptimal grid point over a true positive mode.

The log spacing spans 3-4 decades for typical branch lengths (from $\sim 10^{-4}$ to $\sim 0.5$), so 100 points give roughly 25-33 points per decade. A mode narrower than one grid step would be invisible.

## Impact

Negligible for surfaces with modes wider than one grid step. Unknown for surfaces with sharp peaks. The Dinh and Matsen 2017 K80 counterexample has a broad local max near $t \approx 0.2$ which is well-covered, but non-unimodal GTR models with extreme parameter values could produce sharper modes.

## Proposed action

1. Convergence test: for a fixture surface with a known mode (e.g. a Gaussian bump of width $w$), vary `GRID_SEARCH_POINTS` across $\{20, 50, 100, 200, 500, 1000\}$ and measure the grid's mode-detection accuracy. Establish the minimum grid density that catches a mode of width $w$.
2. Document the result as "grid search catches modes of width $\geq X$ subs/site per decade".
3. If 100 points is insufficient for realistic sharp modes, either raise the constant or switch to a two-pass grid (coarse scan to locate a concave region, then refinement around the best candidate).

## Cross-references

- [packages/treetime/src/commands/optimize/optimize_unified.rs](../../packages/treetime/src/commands/optimize/optimize_unified.rs) (`grid_search_inner`, `grid_search_branch_lengths`, `GRID_SEARCH_POINTS`)
