# Grid search upper bound caps at 0.5 subs/site

> [!WARNING]
> **Needs review.** Action 1 below proposes a rustdoc on `grid_search_branch_lengths`, but the project bans doc comments on items that feed neither CLI help nor a generated schema. The limitation needs another carrier, for example a test that pins the boundary of applicability (action 2) or a KB entry.

> [!IMPORTANT]
> **Decision required.** The proposed actions are alternatives: accept the cap as a known limitation (actions 1 and 2), or extend the grid on the non-unimodal path (action 3). Extending the grid changes optimizer results on saturated edges and needs a v0 comparison before it can be chosen.

## Problem

`fn grid_search_branch_lengths()` [packages/treetime/src/optimize/zero_boundary.rs#L106](../../packages/treetime/src/optimize/zero_boundary.rs#L106) computes the grid upper bound as `max(1.5 * branch_length + one_mutation, GRID_SEARCH_MIN_UPPER)` with `GRID_SEARCH_MIN_UPPER = 0.5` [packages/treetime/src/optimize/zero_boundary.rs#L10](../../packages/treetime/src/optimize/zero_boundary.rs#L10). The floor ensures the grid covers the biologically plausible range even when the current branch length estimate is zero or very small. The upper bound grows proportionally with the input, so for inputs above $0.5 / 1.5 \approx 0.33$ the proportional term dominates.

For a non-unimodal surface whose global maximum lies far beyond $1.5 \cdot \text{input}$ (e.g. a saturated substitution regime with true $t \gg 0.5$), the grid scan cannot reach it. `fn reconcile_zero_boundary()` [packages/treetime/src/optimize/zero_boundary.rs#L40](../../packages/treetime/src/optimize/zero_boundary.rs#L40) would then either return zero (if zero beats every visible grid point) or return the best visible positive mode, both of which are wrong in the saturated regime.

The Brent search interval uses the same bound: `fn brent_bracket()` [packages/treetime/src/optimize/method_brent.rs#L101-L105](../../packages/treetime/src/optimize/method_brent.rs#L101-L105) sets the upper end to `max(1.5 * branch_length + one_mutation, GRID_SEARCH_MIN_UPPER)`, so the Brent methods cannot move a branch beyond this bound in one step either.

## Impact

Negligible for real phylogenetic data. Branch lengths above $0.5$ subs/site are rare (correspond to near-saturation, where the likelihood surface is nearly flat and the optimizer's choice has little downstream impact). Branches above $\approx 1.0$ are typically an artifact of bad data or a misspecified model rather than a true biological signal.

On synthetic fixtures designed to stress the optimizer, and on multi-modal surfaces where the global max is at equilibrium (large $t$), the cap can produce incorrect results without any warning.

## Proposed action

1. Document the cap as an intentional algorithmic limitation with explicit scope: "grid misses modes beyond $t = \max(1.5 \cdot \text{input}, 0.5)$".
2. Add a test with a fixture where the true mode lies beyond the grid range, asserting that the test documents the known boundary of applicability (not that the optimizer finds the mode).
3. Alternative: extend the grid for the non-unimodal path specifically, e.g. cover $[\epsilon, \text{equilibrium\_t}]$ where equilibrium is the $\ln(2) \cdot 10$ cutoff used elsewhere in the code.

## Cross-references

- `fn grid_search_branch_lengths()`, `const GRID_SEARCH_MIN_UPPER` [packages/treetime/src/optimize/zero_boundary.rs](../../packages/treetime/src/optimize/zero_boundary.rs)
- `fn brent_bracket()` [packages/treetime/src/optimize/method_brent.rs](../../packages/treetime/src/optimize/method_brent.rs)
