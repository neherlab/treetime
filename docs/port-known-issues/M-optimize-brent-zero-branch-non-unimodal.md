# Brent paths exclude zero branch length for non-unimodal models

## Problem

All three Brent variants (`brent`, `brent-sqrt`, `brent-log`) use a bracket lower bound of `max(min_branch_length, 1e-12)`, which is strictly positive. For non-unimodal models (K80, HKY85, TN93, general GTR), the `is_zero_branch_optimal` derivative shortcut returns false (it only handles unimodal models). This means edges where the true ML optimum is at zero branch length can be biased positive when using a Brent method with a non-unimodal model.

The Newton path does not have this problem: it can converge to values near zero and the grid search fallback includes `is_zero_better_than_grid_best` which explicitly compares zero against the best positive candidate.

## Impact

Affects Brent methods with non-unimodal GTR models only. With the default JC69/F81 models (unimodal), the derivative shortcut correctly identifies zero-optimal edges before Brent is called. The bias is small (lower bound is `1e-12` or `0.01 * one_mutation`), but the optimizer cannot return exactly zero.

Downstream effect on topology cleanup: edges whose true optimum is $t = 0$ under a multimodal model are forced to a tiny positive value. `find_zero_optimal_internal_edges()` then never collects them for collapse, so the optimize loop's topology cleanup misses these edges. This changes default optimize behavior for non-unimodal models.

## Fix

Add a post-Brent zero-comparison in the dispatch layer (`run_optimize_mixed`), analogous to `is_zero_better_than_grid_best` for the grid search path. When `indel_count == 0` and all sites are valid at zero, compare the Brent result against zero and return the better value.

The fix belongs in the dispatch match arms, not in the individual Brent functions, because the zero-comparison logic is shared across all methods.

## Cross-references

- Bracket computation: `packages/treetime/src/commands/optimize/method_brent.rs` (`brent_bracket`)
- Zero-branch shortcut: `packages/treetime/src/commands/optimize/optimize_unified.rs` (`is_zero_branch_optimal`)
- Grid zero-comparison: `packages/treetime/src/commands/optimize/optimize_unified.rs` (`is_zero_better_than_grid_best`)
- Dispatch: `packages/treetime/src/commands/optimize/optimize_unified.rs` (`run_optimize_mixed`)
- Topology collapse: `packages/treetime/src/commands/optimize/run.rs` (`find_zero_optimal_internal_edges`)
- Unimodal flag: `packages/treetime/src/gtr/gtr.rs` (`unimodal_branch_likelihood`)
