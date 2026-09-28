# TreeRegression.optimal_reroot() turns a single-child root into an extra leaf

## Problem

`TreeRegression.optimal_reroot()` reroots with Biopython's `Tree.root_with_outgroup()`, which turns a single-child root into a leaf (see [kb/reports/reroot-single-child-root.md](../reports/reroot-single-child-root.md#root-cause)). `TreeTime.reroot()` removes an undated single-child root before calling it, so `treetime`, `treetime clock`, and the `TreeTime` API are not affected. Direct use of `TreeRegression` is.

On the tree from the report, direct use gives 101 leaves for 100 samples, with the extra leaf `i373496`.

Direct use also fails with `AttributeError: 'Clade' object has no attribute 'bad_branch'` unless the caller sets `bad_branch` on every node: `_optimal_root_along_branch()` reads it, but `TreeRegression` does not set it.

## Open question

`TreeRegression` knows only the tip values of leaves, so it cannot tell a dated root from an undated one. Removing every single-child root inside `optimal_reroot()` would drop a dated root when `TreeTime` calls it. A fix needs a way for the caller to state which nodes carry data, or a check that the root is undated before `TreeRegression` reroots.
