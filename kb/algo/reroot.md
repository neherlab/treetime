# Rerooting Algorithms

Every scored root search runs the clock root search in [`clock/find_best_root/`](../../packages/treetime/src/clock/find_best_root/): ordinary least squares on root-to-tip divergences, with the exact optimum on every edge ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md)). `clock` and `timetree` reroot with it by least squares (`least-squares`, `clock-filter`) or at a fixed rate of zero (`min-dev`); `optimize --reroot=min-dev` uses the same search at rate zero. [`reroot/`](../../packages/treetime/src/reroot/) holds the topology moves that all callers share and the date-free entry point of `optimize`.

## Root Search

The backward and forward clock passes give, on both sides of every edge, the sums $(n, \sum t, \sum t^2, \sum d, \sum d^2, \sum td)$ of the dated tips that are not outliers (`fn clock_regression()`). With the root at fraction $x$ of an edge of length $\ell$, the tips below the edge are $(1 - x)\ell$ and all others $x\ell$ further from the root (`fn ClockSet::shift()`), so the score is a quadratic in $x$ and the parabola through $x = 0, \tfrac12, 1$ gives its exact minimum (`fn optimize_exact()`). Every edge is searched, the current root is the baseline, and the minimum over all candidates wins. A required positive rate is checked on each edge's optimum.

v1:

- Search: `fn find_best_root()` [packages/treetime/src/clock/find_best_root/find_best_root.rs#L16](../../packages/treetime/src/clock/find_best_root/find_best_root.rs#L16), `fn optimize_exact()` [method_exact.rs#L19](../../packages/treetime/src/clock/find_best_root/method_exact.rs#L19)
- Method dispatch: `fn select_root()` [packages/treetime/src/clock/reroot.rs#L79](../../packages/treetime/src/clock/reroot.rs#L79)
- Topology move: `fn apply_root_at_edge()` [packages/treetime/src/reroot/orchestrate.rs#L100](../../packages/treetime/src/reroot/orchestrate.rs#L100), placement rules in [`reroot/placement.rs`](../../packages/treetime/src/reroot/placement.rs)
- Topology primitives: [`packages/treetime-graph/src/reroot.rs`](../../packages/treetime-graph/src/reroot.rs), [`common_ancestor.rs`](../../packages/treetime-graph/src/common_ancestor.rs)

After the move, the clock reruns the backward pass on the new tree instead of carrying the regression state through the move (`fn reroot_clock_tree()`, [packages/treetime/src/clock/reroot.rs#L20](../../packages/treetime/src/clock/reroot.rs#L20)).

## Minimum-Variance Rooting (min-dev)

`min-dev` ranks candidates with `RootObjective::FixedRate(0.0)`:

$$\text{score} = \frac{S_{dd}\, n - S_d^2}{2n^2}$$

where $n$ is the number of contributing tips, $S_d$ the sum of root-to-tip distances and $S_{dd}$ their sum of squares; the score is half the variance of the root-to-tip distances. Sampling dates do not occur in it, so changing them cannot change the selected root while the contributing tips stay the same.

`optimize` has no dates: `fn reroot_min_dev()` gives every tip the same date, so that every tip contributes, and runs the same search ([packages/treetime/src/reroot/orchestrate.rs#L19](../../packages/treetime/src/reroot/orchestrate.rs#L19)).

TreeTime v0 passes a zero slope for `min_dev`, but its `base_regression()` function still subtracts the estimated-rate regression term when scoring candidates, so its executed behavior does not implement its documented minimum-variance objective.

v0 reference and erratum:

- `min_dev` dispatch with `slope = 0`: [`packages/legacy/treetime/treetime/treetime.py`](../../packages/legacy/treetime/treetime/treetime.py) (`reroot`)
- Fixed-slope scoring defect: [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md)
- Approved v1 behavior: [kb/decisions/clock-min-dev-fixed-zero-rate.md](../decisions/clock-min-dev-fixed-zero-rate.md)

## Tip / MRCA Rooting

`--reroot-tips` roots on the branch leading to the MRCA of the named tips, at its midpoint, without scoring. The MRCA is found by path intersection (`common_ancestor`).
