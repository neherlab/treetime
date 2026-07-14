# Grid search upper bound caps at 0.5 subs/site

## Problem

`grid_search_branch_lengths` at [packages/treetime/src/optimize/dispatch.rs](../../packages/treetime/src/optimize/dispatch.rs) computes the grid upper bound as `max(1.5 * branch_length + one_mutation, GRID_SEARCH_MIN_UPPER)` with `GRID_SEARCH_MIN_UPPER = 0.5`. The floor ensures the grid covers the biologically plausible range even when the current branch length estimate is zero or very small. The upper bound grows proportionally with the input, so for inputs above $0.5 / 1.5 \approx 0.33$ the proportional term dominates.

For a non-unimodal surface whose global maximum lies far beyond $1.5 \cdot \text{input}$ (e.g. a saturated substitution regime with true $t \gg 0.5$), the grid scan cannot reach it. `reconcile_zero_boundary` would then either return zero (if zero beats every visible grid point) or return the best visible positive mode, both of which are wrong in the saturated regime.

## Impact

Negligible for real phylogenetic data. Branch lengths above $0.5$ subs/site are rare (correspond to near-saturation, where the likelihood surface is nearly flat and the optimizer's choice has little downstream impact). Branches above $\approx 1.0$ are typically an artifact of bad data or a misspecified model rather than a true biological signal.

On synthetic fixtures designed to stress the optimizer, and on multi-modal surfaces where the global max is at equilibrium (large $t$), the cap can produce incorrect results without any warning.

## Proposed action

1. Document the cap as an intentional algorithmic limitation in the `grid_search_branch_lengths` rustdoc with explicit scope: "grid misses modes beyond $t = \max(1.5 \cdot \text{input}, 0.5)$".
2. Add a test with a fixture where the true mode lies beyond the grid range, asserting that the test documents the known boundary of applicability (not that the optimizer finds the mode).
3. Alternative: extend the grid for the non-unimodal path specifically, e.g. cover $[\epsilon, \text{equilibrium\_t}]$ where equilibrium is the $\ln(2) \cdot 10$ cutoff used elsewhere in the code.

## Cross-references

- [packages/treetime/src/optimize/dispatch.rs](../../packages/treetime/src/optimize/dispatch.rs) (`grid_search_branch_lengths`, `GRID_SEARCH_MIN_UPPER`)

## Related issues

- Source: [kb/issues/N-optimize-grid-search-upper-bound-capped.md](../issues/N-optimize-grid-search-upper-bound-capped.md) -- delete after full resolution
