# Rerooting Algorithms

Generic, scoring-pluggable root search shared across commands. The search algorithm is decoupled from any particular objective through the `RootStats` trait, so the same machinery serves clock-based rerooting (variance-weighted regression, dated tips) and date-free rerooting (divergence variance).

The clock command currently retains its own copy of the search in `clock/find_best_root/`; its migration contract is tracked in [kb/issues/N-reroot-clock-search-duplicates-generic-module.md](../issues/N-reroot-clock-search-duplicates-generic-module.md).

Design rationale and the full derivation live in the proposals:

- [`kb/proposals/reroot-generic-scoring-architecture.md`](../proposals/reroot-generic-scoring-architecture.md)
- [`kb/proposals/optimize-reroot-support.md`](../proposals/optimize-reroot-support.md)

## Generic Root Search

The objective is abstracted by the `RootStats` trait: sufficient statistics that accumulate per-tip contributions up the tree (`leaf`), push across branches (`propagate`), combine across subtrees (`Add`), recover the complementary message at a node (`Sub`), and reduce to a scalar objective (`score`, minimized).

`EdgeCostFn<S>` evaluates the combined statistics at a split fraction `x` along an edge (`x = 0` at the source/parent, `x = 1` at the target/child), splitting the branch variance linearly. The discrete search optimizes the split on every edge and keeps the global minimum, using the current root's score as the baseline so an already-optimal tree is left unchanged. Optimizing each edge over `[0, 1]` covers rooting at any existing node as a split endpoint.

Type bounds shrink from the clock-specific `ClockNode`/`ClockEdge` to `GraphNode`/`GraphEdge`; statistics are passed in a `BTreeMap<GraphEdgeKey, (S, S)>` rather than read from node/edge payloads, so the search works on payloads (e.g. the optimize graph) that carry no message fields.

v1:

- Trait: [`packages/treetime/src/reroot/traits.rs`](../../packages/treetime/src/reroot/traits.rs)
- Cost function: [`packages/treetime/src/reroot/cost_function.rs`](../../packages/treetime/src/reroot/cost_function.rs)
- Search: [`packages/treetime/src/reroot/search.rs`](../../packages/treetime/src/reroot/search.rs), [`split.rs`](../../packages/treetime/src/reroot/split.rs)
- Orchestration (topology mutation + domain fixup callback): [`packages/treetime/src/reroot/orchestrate.rs`](../../packages/treetime/src/reroot/orchestrate.rs)
- Topology primitives: [`packages/treetime-graph/src/reroot.rs`](../../packages/treetime-graph/src/reroot.rs), [`common_ancestor.rs`](../../packages/treetime-graph/src/common_ancestor.rs)

## Divergence Rooting (min-dev)

`DivStats` scores a root position by the weighted variance of root-to-tip distances:

$$\text{score} = \frac{S_{dd}\, n - S_d^2}{2n}$$

where $n$ is the summed inverse branch variance over tips, $S_d$ the variance-weighted sum of root-to-tip distances, and $S_{dd}$ the weighted sum of squares. This requires no dates: every tip contributes regardless of whether it carries one.

These are the divergence components of the clock sufficient statistics (`ClockSet`) with the date-gated `norm` replaced by an always-populated `count`; `propagate` is `ClockSet::propagate_averages` with all time terms identically zero. The objective is exactly v0 `min_dev` rooting: least-squares regression with a fixed zero slope, whose chi-squared with no dates reduces to $\tfrac12(S_{dd} - S_d^2/n)$.

Variance model (`VarianceModel`, defaults `factor = 0`, `offset = 0`, `offset_leaf = 1`): internal branches have variance 0, leaf branches variance 1. This reproduces v0's no-covariation `TreeRegression` weighting, under which the score is the plain variance of root-to-tip distances.

v0 reference:

- `min_dev` dispatch with `slope = 0`: [`packages/legacy/treetime/treetime/treetime.py`](../../packages/legacy/treetime/treetime/treetime.py) (`reroot`)
- Fixed-slope chi-squared: [`packages/legacy/treetime/treetime/treeregression.py`](../../packages/legacy/treetime/treetime/treeregression.py) (`base_regression`)
- No-covariation variance: [`packages/legacy/treetime/treetime/clock_tree.py`](../../packages/legacy/treetime/treetime/clock_tree.py) (`setup_TreeRegression`)

v1:

- Statistics: [`packages/treetime/src/reroot/div_stats.rs`](../../packages/treetime/src/reroot/div_stats.rs)
- Two-pass traversal: [`packages/treetime/src/reroot/div_stats_traversal.rs`](../../packages/treetime/src/reroot/div_stats_traversal.rs)

## Tip / MRCA Rooting

`--reroot-tips` roots on the branch leading to the MRCA of the named tips, at its midpoint, without scoring. The MRCA is found by path intersection (`common_ancestor`).
