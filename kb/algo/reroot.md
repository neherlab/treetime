# Rerooting Algorithms

Generic, scoring-pluggable root search currently used for date-free divergence rooting. The search algorithm is decoupled from a particular objective through the `RootStats` trait and is designed to accept clock sufficient statistics as well.

Clock rerooting still uses a parallel search in `clock/find_best_root/`; it does not yet share the generic engine. Its migration contract is tracked in [kb/issues/N-reroot-clock-search-duplicates-generic-module.md](../issues/N-reroot-clock-search-duplicates-generic-module.md).

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

`DivStats` scores a root position by half the weighted residual sum around the weighted mean root-to-tip distance:

$$\text{score} = \frac{S_{dd}\, n - S_d^2}{2n}$$

where $n$ is the summed inverse branch variance over tips, $S_d$ the variance-weighted sum of root-to-tip distances, and $S_{dd}$ the weighted sum of squares. Dividing this score by $n/2$ gives the weighted variance. This requires no dates: every tip contributes regardless of whether it carries one.

These are the divergence components of the clock sufficient statistics (`ClockSet`) with the date-gated `norm` replaced by an always-populated `count`; `propagate` is `ClockSet::propagate_averages` with all time terms identically zero. The objective implements the documented `min_dev` meaning: least-squares regression with a fixed zero slope, whose unnormalized chi-squared reduces to $\tfrac12(S_{dd} - S_d^2/n)$.

TreeTime v0 passes a zero slope for `min_dev`, but its `base_regression()` function still subtracts the estimated-rate regression term when scoring candidates. Its executed behavior therefore does not implement its documented minimum-variance objective. See [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md).

Variance model (`VarianceModel`, defaults `factor = 0`, `offset = 0`, `offset_leaf = 1`): internal branches have variance 0, leaf branches variance 1. This reproduces v0's no-covariation `TreeRegression` weighting, under which the score is proportional to the plain variance of root-to-tip distances.

## Clock Rooting (min-dev)

Clock `min-dev` uses the clock sufficient statistics and ranks candidates with `RootObjective::FixedRate(0.0)`:

$$\text{score}_{\mathrm{clock}} = \frac{S_{dd}\, n - S_d^2}{2n^2}$$

where $n$, $S_d$, and $S_{dd}$ have the same meanings as above. Sampling-date values do not occur in this expression, so changing them cannot change the selected root while the contributing tip set and variance model remain fixed.

The generic and clock scores differ by a factor of $n$. Under the default no-covariation model, $n$ is the constant number of contributing tips, so both choose the same root. Under a covariation model, effective precision can vary with root position; clock `min-dev` retains the normalized WLS objective while generic `DivStats` retains its date-free weighted residual-sum objective.

v0 reference and erratum:

- `min_dev` dispatch with `slope = 0`: [`packages/legacy/treetime/treetime/treetime.py`](../../packages/legacy/treetime/treetime/treetime.py) (`reroot`)
- Fixed-slope scoring defect: [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md)
- Approved v1 behavior: [kb/decisions/clock-min-dev-fixed-zero-rate.md](../decisions/clock-min-dev-fixed-zero-rate.md)
- No-covariation variance: [`packages/legacy/treetime/treetime/clock_tree.py`](../../packages/legacy/treetime/treetime/clock_tree.py) (`setup_TreeRegression`)

v1:

- Statistics: [`packages/treetime/src/reroot/div_stats.rs`](../../packages/treetime/src/reroot/div_stats.rs)
- Two-pass traversal: [`packages/treetime/src/reroot/div_stats_traversal.rs`](../../packages/treetime/src/reroot/div_stats_traversal.rs)
- Clock dispatch: [`packages/treetime/src/clock/reroot.rs`](../../packages/treetime/src/clock/reroot.rs)

## Tip / MRCA Rooting

`--reroot-tips` roots on the branch leading to the MRCA of the named tips, at its midpoint, without scoring. The MRCA is found by path intersection (`common_ancestor`).
