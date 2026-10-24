# Rerooting Algorithms

The top-level module [`reroot/`](../../packages/treetime/src/reroot/) holds every rooting method of `clock`, `timetree` and `optimize`. Choosing a root and moving it are separate steps:

- **Root search**: `fn find_root()` [packages/treetime/src/reroot/find_root.rs#L15](../../packages/treetime/src/reroot/find_root.rs#L15) reads the tree, and the dates for the date-based methods, and returns a root position: an edge and a split $x \in [0, 1]$ measured from the parent end of the edge, or "keep the current root". It never changes the tree. One function per method:
  - `least-squares` and `min-dev`: the edge search, `fn search_edges()` [packages/treetime/src/reroot/edge_search.rs#L19](../../packages/treetime/src/reroot/edge_search.rs#L19)
  - `oldest` and `--outgroup`: the named-node target, `fn named_node_point()` [packages/treetime/src/reroot/named_node.rs#L13](../../packages/treetime/src/reroot/named_node.rs#L13)
  - `midpoint`: `fn midpoint_root()` [packages/treetime/src/reroot/midpoint.rs#L14](../../packages/treetime/src/reroot/midpoint.rs#L14)
- **Topology change**: `fn apply_root_at_edge()` [packages/treetime/src/reroot/orchestrate.rs#L20](../../packages/treetime/src/reroot/orchestrate.rs#L20) turns the root position into a node or an edge split (`fn RootTarget::on_edge()` [packages/treetime/src/reroot/placement.rs#L62](../../packages/treetime/src/reroot/placement.rs#L62)) and moves the root. The topology primitives are in [`packages/treetime-graph/src/reroot.rs`](../../packages/treetime-graph/src/reroot.rs), and `struct SplitFraction` ([`packages/treetime-graph/src/split_fraction.rs`](../../packages/treetime-graph/src/split_fraction.rs)) keeps every split finite and inside $[0, 1]$
- **Driver**: `fn reroot_at_best_root()` [packages/treetime/src/reroot/driver.rs#L16](../../packages/treetime/src/reroot/driver.rs#L16) searches, removes an undated stem root (a root with one child) and searches again when the root moves, applies the move, and gives the new nodes and edges default clock inputs. A dated stem root stays and becomes a dated tip

The clock fit (`fn fit_clock_model()`, [`clock/clock_regression.rs`](../../packages/treetime/src/clock/clock_regression.rs)) only fits. Every caller reroots first and then fits. The clock regression files never import from `reroot/`; `reroot/` reads the regression sums from `clock`.

## Date-Based and Date-Free Methods

`fn RerootSpec::uses_dates()` [packages/treetime/src/reroot/params.rs#L15](../../packages/treetime/src/reroot/params.rs#L15) splits the methods:

- **Date-based**: `least-squares` and `oldest`. `clock` reroots before and after its outlier filter, and `timetree` reroots before its first clock fit and twice in its pre-loop ([kb/issues/M-timetree-reroot-schedule-diverges-from-v0.md](../issues/M-timetree-reroot-schedule-diverges-from-v0.md))
- **Date-free**: `min-dev`, `--outgroup` and `midpoint`. They root once, on the starting branch lengths, through `fn reroot_date_free()` [packages/treetime/src/reroot/date_free.rs#L15](../../packages/treetime/src/reroot/date_free.rs#L15): in `clock` and `timetree` before the first clock fit, on the input lengths, and in `optimize` right after the initial guess. The root then stays, and every later clock fit of the run still happens. In `timetree`, a clock refit on the kept root that gives no positive regression rate starts from the rate of the Gaussian timetree, as the first clock fit and the clock updates of the refinement loop do (`fn gaussian_start_rate()`). The search counts every tip, dated or not, and leaves no clock-filter outlier out. An `is_dated` predicate drives the stem-root check and the check that rejects an undated inner node that the move turns into a leaf: `clock` and `timetree` pass the run's dates, `optimize` passes "no node has a date". `min-dev` and `midpoint` place the root from branch lengths and stop with an input error when an edge has a missing, negative or non-finite length; `--outgroup` needs no lengths ([kb/decisions/reroot-date-free-methods-root-once.md](../decisions/reroot-date-free-methods-root-once.md))

## Edge Search

The backward and forward clock passes give, on both sides of every edge, the sums $(n, \sum t, \sum t^2, \sum d, \sum d^2, \sum td)$ of the tips that contribute (`fn clock_regression()`): the dated tips that are not outliers for `least-squares`, every tip for `min-dev`. With the root at fraction $x$ of an edge of length $\ell$, the tips below the edge are $(1 - x)\ell$ and all others $x\ell$ further from the root (`fn ClockSet::shift()`), so the score is a quadratic in $x$ and the parabola through $x = 0, \tfrac12, 1$ gives its exact minimum (`fn optimize_exact()`). Every edge is searched, the current root is the baseline, and the minimum over all candidates wins. A required positive rate is checked on each edge's optimum. The search is ordinary least squares on root-to-tip divergences ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md)).

## Root Placement

`fn RootTarget::on_edge()` turns a root position into the new root. It has two placement modes, chosen by the method (`fn RerootSpec::placement()`):

- **Edge-end band** (`least-squares`, `min-dev`): a split $x < 10^{-5}$ roots at the parent node, and $x > 1 - 10^{-5}$ at the child node
- **Exact** (`oldest`, `--outgroup`, `midpoint`): positions computed from branch lengths move to an edge end only when their distance from that end, $x\ell$ or $(1 - x)\ell$, is zero within 4 ulps of $\ell$ (`fn coordinates_coincide()`). Biopython, and so augur, also place these positions exactly

In both modes, any other $x$ adds a new node that splits the edge at $x$, with the lengths $x\ell$ above and $(1 - x)\ell$ below it. Without edge splitting (`split_edge = false`), the root snaps to the nearer end: the parent for $x < \tfrac12$, the child otherwise.

**The edge search places a root within $10^{-5}$ of an edge end at the end node.** The score along an edge is the quadratic $q(x) = a x^2 + b x + c$, so a root at the end node instead of at an optimum $x^* < 10^{-5}$ raises the score by $a x^{*2} < 10^{-10} a$. A new node at $x^*$ would add a root branch of length $x^* \ell$ that changes no result. The band is measured on the split fraction, the same band as v0 (below). The search returns exactly 0 or 1 when the optimum lies beyond an edge end, so the band only decides optima that lie inside the edge.

The band also keeps the root when a later search finds the same root again. `clock` with the clock filter searches once before and once after the filter. When the second search finds the current root, rounding puts its optimum a small fraction inside one of the root's edges, larger than an ulp (below $5 \cdot 10^{-7}$ on `data/dengue/20`, `data/flu/h3n2/20` and `data/zika/86`). A split there would add a new root node, rename the root, and change branch lengths in the last digits when the old root is merged away.

**A leaf never becomes the root.** On a terminal edge the exact optimum often lies past the tip: the vertex of the parabola is beyond $x = 1$, so the clamp gives $x = 1$. With edge splitting, the reroot then inserts a new root at the tip end, with a branch of length 0 to the tip and the full length $\ell$ on the other side. Without edge splitting, the parent end becomes the root. A tip as root would give the root the tip's exact date while the tips of its subtree are older, so the time inference would find an empty time distribution at the root and leave the tip undated. `fn apply_root_at_edge()` rejects a leaf as the new root with an internal error.

The band and the leaf rule give the v0 result. v0 measures the split from the node toward its parent, so its $x$ is $1 - x$ of v1:

- v0 scores the node end of a terminal edge as infinite, so a tip is never a candidate root node ([packages/legacy/treetime/treetime/treeregression.py#L397](../../packages/legacy/treetime/treetime/treeregression.py#L397))
- A split with $x < 10^{-5}$ roots at the edge's node through `root_with_outgroup(node)`, and a split with $x > 1 - 10^{-5}$ at its parent ([packages/legacy/treetime/treetime/treeregression.py#L443-L463](../../packages/legacy/treetime/treetime/treeregression.py#L443-L463))
- For a terminal outgroup, Biopython's `root_with_outgroup()` creates a new bifurcating root with a 0-length branch to the outgroup ([Bio/Phylo/BaseTree.py#L843-L846](https://github.com/biopython/biopython/blob/d7e4b8b19399668b09442a5b35765d9186b5f665/Bio/Phylo/BaseTree.py#L843-L846), Biopython 1.88)

## Minimum-Variance Rooting (min-dev)

`min-dev` ranks candidates with `RootObjective::FixedRate(0.0)`, whose score is `fn ClockSet::chisq_fixed_rate()` at rate zero:

$$\text{score} = \frac{1}{2}\left(S_{dd} - \frac{S_d^2}{n}\right)$$

where $n$ is the number of contributing tips, $S_d$ the sum of root-to-tip distances and $S_{dd}$ their sum of squares; the score is half the sum of squared deviations of the root-to-tip distances from their mean. Sampling dates do not occur in it, and the date-free search counts every tip, so the dates cannot change the selected root.

TreeTime v0 passes a zero slope for `min_dev`, but its `base_regression()` function still subtracts the estimated-rate regression term when scoring candidates, so its executed behavior does not implement its documented minimum-variance objective.

v0 reference and erratum:

- `min_dev` dispatch with `slope = 0`: [`packages/legacy/treetime/treetime/treetime.py`](../../packages/legacy/treetime/treetime/treetime.py) (`reroot`)
- Fixed-slope scoring defect: [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md)
- Approved v1 behavior: [kb/decisions/clock-min-dev-fixed-zero-rate.md](../decisions/clock-min-dev-fixed-zero-rate.md)

## Rooting on a Named Node (oldest, --outgroup)

`--outgroup` names tips, and the named node $M$ is their most recent common ancestor, found by path intersection (`common_ancestor`). `--reroot=oldest` names the oldest dated tip. The root goes to the middle of the branch above $M$ in the unrooted tree, without scoring. With $R$ the parent of $M$:

- $R$ is not the root, or $R$ does not have exactly two children: the edge `R -> M` at the split $\tfrac12$
- $R$ is the root with the two children $M$ (branch length $a$) and $S$ (branch length $b$): the branch of $M$ in the unrooted tree is `R -> M` plus `R -> S`, of length $a + b$. The root lies $(a+b)/2$ from $M$: on `R -> M` at $x = (a-b)/(2a)$ when $a > b$, on `R -> S` at $x = (b-a)/(2b)$ when $b > a$, and at $R$ when $a = b$ or a length is missing or invalid

So a second reroot on the same node changes nothing, also when a user reroots the output of an earlier run. v0 halves the node's current branch on every reroot ([kb/v0-errata/reroot-named-node-halves-branch-on-each-call.md](../v0-errata/reroot-named-node-halves-branch-on-each-call.md), [kb/decisions/reroot-named-node-at-unrooted-branch-middle.md](../decisions/reroot-named-node-at-unrooted-branch-middle.md)).

## Midpoint Rooting

`midpoint` roots at the middle of the longest tip-to-tip path, using only branch lengths. This is the behavior of `augur refine --root mid_point`: augur calls Biopython's `Tree.root_at_midpoint()` before TreeTime runs, then keeps that root ([augur/refine.py#L134-L137](https://github.com/nextstrain/augur/blob/58f329042e7fc631ca275cf2df44c8720f80949f/augur/refine.py#L134-L137), [#L285-L287](https://github.com/nextstrain/augur/blob/58f329042e7fc631ca275cf2df44c8720f80949f/augur/refine.py#L285-L287)). Biopython reroots at every tip to find the two farthest tips, so its time grows with the square of the tip count ([Bio/Phylo/BaseTree.py#L904-L938](https://github.com/biopython/biopython/blob/d7e4b8b19399668b09442a5b35765d9186b5f665/Bio/Phylo/BaseTree.py#L904-L938), Biopython 1.88). `fn midpoint_root()` is linear in the number of nodes.

The tree is read as unrooted: its tips are the leaves, plus the root when it has exactly one child. With $\ell(e)$ the length of edge $e$:

- A postorder pass gives every node $v$ its height $h(v)$, the longest distance down to a leaf, and the leaf that reaches it. A leaf has $h = 0$, and an inner node $h(v) = \max_{e = (v \to c)} \big(h(c) + \ell(e)\big)$
- A node with at least two children closes the longest path whose top node is $v$, of length $h_1(v) + h_2(v)$, the two largest child values. A node with one child only passes its height up, except the root: a root with one child closes the path from itself down to its deepest leaf, of length $h(\text{root})$
- $D$ is the longest such path, with end tips $u$ (on the longer arm) and $w$. For ties of $D$, the first maximum in postorder and in the order of the children wins, so the result depends only on the input file
- From $u$, the walk up the path covers the distance $s$. The root lies on the first edge $e = (p \to c)$ with $\ell(e) > 0$ and $s + \ell(e) \ge D/2$, at $x = (s + \ell(e) - D/2) / \ell(e)$ from the parent end. When the two arms are equal, this gives $x = 0$ on the last edge of the arm, so the root is the top node

The position is placed exactly. An undated stem root is a tip of the unrooted tree, so the midpoint never lies on it and the root moves; the driver then removes the stem and searches again on the tree without it. A dated stem root stays and becomes a tip, and the search has already counted it as one. With fewer than two tips, or $D = 0$, the root stays. Dates and outliers are not read.

Midpoint rooting and `min-dev` both assume a strict clock with all tips sampled at the same time. Midpoint rooting minimizes the largest root-to-tip distance and uses only the two farthest tips; `min-dev` uses all tips, and was more accurate than midpoint rooting in the simulations of [Mai, Sayyari and Mirarab 2017](https://doi.org/10.1371/journal.pone.0182238).

The golden master `test_gm_reroot_midpoint` compares v1 with Biopython 1.88 `Tree.root_at_midpoint()` on five datasets ([`reroot/__tests__/__fixtures__/gm_reroot_midpoint_capture`](../../packages/treetime/src/reroot/__tests__/__fixtures__/gm_reroot_midpoint_capture)).
