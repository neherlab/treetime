# Rerooting Algorithms

Every scored root search runs the clock root search in [`clock/find_best_root/`](../../packages/treetime/src/clock/find_best_root/): ordinary least squares on root-to-tip divergences, with the exact optimum on every edge ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md)). `clock` and `timetree` reroot with it by least squares (`least-squares`, `clock-filter`) or at a fixed rate of zero (`min-dev`); `optimize --reroot=min-dev` uses the same search at rate zero. [`reroot/`](../../packages/treetime/src/reroot/) holds the topology move that the clock driver uses; `optimize` reroots through the same driver.

## Root Search

The backward and forward clock passes give, on both sides of every edge, the sums $(n, \sum t, \sum t^2, \sum d, \sum d^2, \sum td)$ of the dated tips that are not outliers (`fn clock_regression()`). With the root at fraction $x$ of an edge of length $\ell$, the tips below the edge are $(1 - x)\ell$ and all others $x\ell$ further from the root (`fn ClockSet::shift()`), so the score is a quadratic in $x$ and the parabola through $x = 0, \tfrac12, 1$ gives its exact minimum (`fn optimize_exact()`). Every edge is searched, the current root is the baseline, and the minimum over all candidates wins. A required positive rate is checked on each edge's optimum.

v1:

- Search: `fn find_best_root()` [packages/treetime/src/clock/find_best_root/find_best_root.rs#L16](../../packages/treetime/src/clock/find_best_root/find_best_root.rs#L16), `fn optimize_exact()` [method_exact.rs#L8](../../packages/treetime/src/clock/find_best_root/method_exact.rs#L8)
- Method dispatch: `fn select_root()` [packages/treetime/src/clock/reroot.rs#L111](../../packages/treetime/src/clock/reroot.rs#L111)
- Topology move: `fn apply_root_at_edge()` [packages/treetime/src/reroot/orchestrate.rs#L19](../../packages/treetime/src/reroot/orchestrate.rs#L19), placement rules in [`reroot/placement.rs`](../../packages/treetime/src/reroot/placement.rs)
- Topology primitives: [`packages/treetime-graph/src/reroot.rs`](../../packages/treetime-graph/src/reroot.rs), [`common_ancestor.rs`](../../packages/treetime-graph/src/common_ancestor.rs)

One driver serves every caller (`fn reroot_at_best_root()`, [packages/treetime/src/clock/reroot.rs#L20](../../packages/treetime/src/clock/reroot.rs#L20)): it searches, removes an undated stem and searches again when the root moves, applies the move and gives the new nodes and edges default clock inputs. The clock then reruns the backward pass on the new tree instead of carrying the regression state through the move.

## Root Placement

The search returns an edge and a split $x \in [0, 1]$, measured from the parent end ($x = 0$) to the child end ($x = 1$) of the edge. `fn RootTarget::on_edge()` turns it into the new root:

- $x < 10^{-5}$: the parent node becomes the root
- $x > 1 - 10^{-5}$: the child node becomes the root, unless the child is a leaf (below)
- Other $x$: a new node splits the edge at $x$, with the lengths $x\ell$ above and $(1 - x)\ell$ below it. Without edge splitting (`split_edge = false`), the root snaps to the nearer end: the parent for $x < \tfrac12$, the child otherwise

**A root within $10^{-5}$ of an edge end goes to the end node.** The score along an edge is the quadratic $q(x) = a x^2 + b x + c$, so a root at the end node instead of at an optimum $x^* < 10^{-5}$ raises the score by $a x^{*2} < 10^{-10} a$. A new node at $x^*$ would add a root branch of length $x^* \ell$ that changes no result. The band is measured on the split fraction, the same band as v0 (below). The search returns exactly 0 or 1 when the optimum lies beyond an edge end, so the band only decides optima that lie inside the edge.

The band also keeps the root when a later search finds the same root again. `clock` with the clock filter searches once before and once after the filter. When the second search finds the current root, rounding puts its optimum a small fraction inside one of the root's edges, larger than an ulp (below $5 \cdot 10^{-7}$ on `data/dengue/20`, `data/flu/h3n2/20` and `data/zika/86`). A split there would add a new root node, rename the root, and change branch lengths in the last digits when the old root is merged away.

**A leaf never becomes the root.** On a terminal edge the exact optimum often lies past the tip: the vertex of the parabola is beyond $x = 1$, so the clamp gives $x = 1$. With edge splitting, the reroot then inserts a new root at the tip end, with a branch of length 0 to the tip and the full length $\ell$ on the other side. Without edge splitting, the parent end becomes the root. A tip as root would give the root the tip's exact date while the tips of its subtree are older, so the time inference would find an empty time distribution at the root and leave the tip undated. `fn apply_root_at_edge()` rejects a leaf as the new root with an internal error.

Both rules give the v0 result. v0 measures the split from the node toward its parent, so its $x$ is $1 - x$ of v1:

- v0 scores the node end of a terminal edge as infinite, so a tip is never a candidate root node ([packages/legacy/treetime/treetime/treeregression.py#L397](../../packages/legacy/treetime/treetime/treeregression.py#L397))
- A split with $x < 10^{-5}$ roots at the edge's node through `root_with_outgroup(node)`, and a split with $x > 1 - 10^{-5}$ at its parent ([packages/legacy/treetime/treetime/treeregression.py#L443-L463](../../packages/legacy/treetime/treetime/treeregression.py#L443-L463))
- For a terminal outgroup, Biopython's `root_with_outgroup()` creates a new bifurcating root with a 0-length branch to the outgroup ([Bio/Phylo/BaseTree.py#L843-L846](https://github.com/biopython/biopython/blob/d7e4b8b19399668b09442a5b35765d9186b5f665/Bio/Phylo/BaseTree.py#L843-L846), Biopython 1.88)

v1: `fn RootTarget::on_edge()` [packages/treetime/src/reroot/placement.rs#L30](../../packages/treetime/src/reroot/placement.rs#L30)

## Minimum-Variance Rooting (min-dev)

`min-dev` ranks candidates with `RootObjective::FixedRate(0.0)`:

$$\text{score} = \frac{S_{dd}\, n - S_d^2}{2n^2}$$

where $n$ is the number of contributing tips, $S_d$ the sum of root-to-tip distances and $S_{dd}$ their sum of squares; the score is half the variance of the root-to-tip distances. Sampling dates do not occur in it, so changing them cannot change the selected root while the contributing tips stay the same.

`optimize` has no dates: `fn reroot_undated()` gives every tip the same date, so that every tip contributes, and runs the same driver at rate zero for `min-dev` and `--reroot-tips` ([packages/treetime/src/clock/reroot.rs#L53](../../packages/treetime/src/clock/reroot.rs#L53)).

TreeTime v0 passes a zero slope for `min_dev`, but its `base_regression()` function still subtracts the estimated-rate regression term when scoring candidates, so its executed behavior does not implement its documented minimum-variance objective.

v0 reference and erratum:

- `min_dev` dispatch with `slope = 0`: [`packages/legacy/treetime/treetime/treetime.py`](../../packages/legacy/treetime/treetime/treetime.py) (`reroot`)
- Fixed-slope scoring defect: [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md)
- Approved v1 behavior: [kb/decisions/clock-min-dev-fixed-zero-rate.md](../decisions/clock-min-dev-fixed-zero-rate.md)

## Tip / MRCA Rooting

`--reroot-tips` roots on the branch leading to the MRCA of the named tips, at its midpoint, without scoring. The MRCA is found by path intersection (`common_ancestor`).
