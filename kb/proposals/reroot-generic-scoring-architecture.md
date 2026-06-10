# Generic root search architecture with pluggable scoring

The root search algorithm in `clock/find_best_root/` is fused with clock regression scoring (`ClockSet.chisq()`). Two independent problems prevent reuse for divergence-only rooting in optimize:

1. Type system: `find_best_root` requires `ClockNode` + `ClockEdge` trait bounds. Optimize uses `GraphAncestral` with `NodeAncestral`/`EdgeAncestral` (`payload/ancestral.rs`), which do not implement these traits. The existing reroot functions cannot be called from optimize at all -- it is a compile error.
2. Mathematical: even if the trait bounds were satisfied, `ClockSet.chisq()` produces NaN when no tips have dates (see proof below).

This proposal extracts the generic search algorithm and makes the scoring function pluggable via a `RootStats` trait, following the `argmin::CostFunction` pattern already used for 1D optimization. The generic search requires only `GraphNode` + `GraphEdge` bounds, solving problem 1. The pluggable `DivStats` scoring uses a formula that works without dates, solving problem 2.

## Problem: ClockSet breaks without dates

`ClockSet` (`payload/clock_set.rs`) aggregates six sufficient statistics for weighted-least-squares regression of divergence against time:

| Field     | Meaning                                 | Value without dates |
| --------- | --------------------------------------- | ------------------- |
| `t_sum`   | $\sum t_i / \sigma_i^2$                 | 0                   |
| `tsq_sum` | $\sum t_i^2 / \sigma_i^2$               | 0                   |
| `d_sum`   | $\sum d_i / \sigma_i^2$                 | nonzero             |
| `dsq_sum` | $\sum d_i^2 / \sigma_i^2$               | nonzero             |
| `dt_sum`  | $\sum d_i t_i / \sigma_i^2$             | 0                   |
| `norm`    | $\sum 1 / \sigma_i^2$ (dated tips only) | **0**               |

The `norm` field is gated on `date.is_some()` (`clock_set.rs:26,49`). When no tips have dates, `norm = 0` everywhere. The scoring function `chisq()` (`clock_set.rs:118-125`) computes:

$$\chi^2 = \frac{1}{2n}\left(S_{dd}\,n - S_d^2 - \frac{(S_{dt}\,n - S_d\,S_t)^2}{\Delta}\right)$$

where $n$ = `norm`, $\Delta$ = `determinant()` = $S_{tt}\,n - S_t^2$. Both $n$ and $\Delta$ are zero, producing $0/0$ = NaN.

This is not a fixable edge case. The formula regresses divergence against time. Without time data, the regression is undefined. A divergence-only objective requires a different formula over different sufficient statistics.

## Algorithm decomposition

The root-finding algorithm has four layers. Only layer 3 changes between objectives.

### Layer 1: topology (already generic)

`treetime-graph/src/reroot.rs`. Domain-agnostic graph operations:

- `split_edge`: insert a new node at position x along an edge
- `apply_reroot_topology`: invert edges on the path from old root to new root
- `remove_node_if_trivial`: merge a degree-2 node with its neighbors

No changes needed. Already parameterized over `N: GraphNode, E: GraphEdge`.

### Layer 2: search (to be made generic)

Currently in `clock/find_best_root/`. Two phases:

1. Discrete scan: iterate every node, evaluate score, keep the best. `find_best_root.rs`.
2. Continuous refinement: on the best node's parent edge, optimize split position $x \in [0, 1]$ via 1D minimizer. `find_best_split.rs` dispatches to `method_brent.rs`, `method_golden_section.rs`, or `method_grid_search.rs`.

The search logic does not depend on what the score means. It calls `score()` and picks the minimum. The 1D methods use `argmin::CostFunction` -- already trait-based.

Currently hardwired to `ClockNode`/`ClockEdge` bounds because `BranchPointCostFunction` reads `ClockSet` fields directly. Generalizing the cost function generalizes the search.

### Layer 3: scoring (to be made pluggable)

Currently in `clock/find_best_root/cost_function.rs`. `BranchPointCostFunction` holds `to_parent: ClockSet` and `to_child: ClockSet` (per-edge sufficient statistics computed during clock regression), and evaluates the combined score at split position x:

```
child_stats = to_parent.propagate(bl * (1 - x), var * (1 - x))
parent_stats = to_child.propagate(bl * x, var * x)
combined = parent_stats + child_stats
score = combined.chisq()
```

The structure (propagate, combine, score) is generic. The statistics type and scoring formula are specific to clock regression. A `RootStats` trait captures this pattern.

### Layer 4: domain fixup (varies per consumer)

After topology changes, domain-specific edge data may need updating. `clock/reroot.rs:314-339` swaps `to_parent`/`to_child` ClockSet messages on inverted edges. For divergence-only reroot, no persistent edge data exists -- fixup is a no-op.

## RootStats trait

The trait captures the sufficient statistics pattern shared by all root-scoring objectives:

```rust
pub trait RootStats: Clone + Default + Add<Output = Self> + Send + Sync {
    fn leaf(time: Option<f64>, branch_length: f64, variance: f64) -> Self;
    fn propagate(&self, branch_length: f64, variance: f64) -> Self;
    fn score(&self) -> f64;
}
```

- `leaf`: what a tip contributes toward its parent. The `time` parameter is `Some(date)` for dated tips, `None` otherwise. Implementations that don't use dates ignore it.
- `propagate`: push statistics through a branch of given length and variance. Returns the accumulated statistics as seen from the other end of the branch. Matches `ClockSet::propagate_averages`.
- `score`: scalar objective value. Lower is better. The search minimizes this.
- `Add`: combine statistics from independent subtrees (left child + right child at an internal node).

### ClockStats (existing ClockSet, implementing RootStats)

Wraps the existing `ClockSet` with `RootStats` implementation:

- `leaf`: existing `leaf_contribution_to_parent`, gates `norm` on `time.is_some()`
- `propagate`: existing `propagate_averages`
- `score`: existing `chisq()` (WLS regression residual)

No behavioral change for clock and timetree commands.

### DivStats (new, divergence-only)

Sufficient statistics for root-to-tip distance variance:

```rust
pub struct DivStats {
    count: f64,
    d_sum: f64,
    dsq_sum: f64,
}
```

These are the divergence components of `ClockSet` with `count` replacing `norm`. The difference: `count` is always incremented for every tip, regardless of whether it has a date.

Implementation:

- `leaf`: `count = 1/var`, `d_sum = bl/var`, `dsq_sum = bl^2/var`. Always contributes (no date gate).
- `propagate`: same formulas as `ClockSet::propagate_averages` with time terms zeroed out:

```rust
fn propagate(&self, bl: f64, var: f64) -> Self {
    let denom = 1.0 / (1.0 + var * self.count);
    DivStats {
        count: self.count * denom,
        d_sum: (self.d_sum + bl * self.count) * denom,
        dsq_sum: self.dsq_sum
            + 2.0 * bl * self.d_sum
            + bl.powi(2) * self.count
            - var * (self.d_sum + bl * self.count).powi(2) * denom,
    }
}
```

- `score`: weighted variance of root-to-tip distances:

$$\text{score} = \frac{S_{dd} \cdot n - S_d^2}{2n}$$

where $n$ = `count`, $S_d$ = `d_sum`, $S_{dd}$ = `dsq_sum`. This is `ClockSet.chisq()` without the time regression term (the $(S_{dt}\,n - S_d\,S_t)^2 / \Delta$ subtrahend disappears when all time terms are zero). `count` is always positive, so no division by zero.

## Generic cost function

`EdgeCostFn<S: RootStats>` replaces `BranchPointCostFunction`:

```rust
pub struct EdgeCostFn<S: RootStats> {
    to_parent: S,
    to_child: S,
    branch_length: f64,
    branch_variance: f64,
    is_leaf: bool,
    leaf_time: Option<f64>,
}
```

Implements `argmin::CostFunction`:

```rust
impl<S: RootStats> CostFunction for &EdgeCostFn<S> {
    type Param = f64;
    type Output = f64;

    fn cost(&self, x: &f64) -> Result<f64, Error> {
        if *x < 0.0 || *x > 1.0 {
            return Ok(f64::INFINITY);
        }
        let child_contrib = if self.is_leaf {
            S::leaf(self.leaf_time, self.branch_length * (1.0 - x), ...)
        } else {
            self.to_parent.propagate(self.branch_length * (1.0 - x), ...)
        };
        let parent_contrib = self.to_child.propagate(self.branch_length * x, ...);
        Ok((parent_contrib + child_contrib).score())
    }
}
```

The 1D methods (`method_brent.rs`, `method_golden_section.rs`, `method_grid_search.rs`) take `impl CostFunction<Param = f64, Output = f64>`. They already work with any cost function. No changes needed.

## Generic search functions

```rust
pub fn find_best_root<S: RootStats>(
    graph: &Graph<N, E, D>,
    edge_stats: &BTreeMap<GraphEdgeKey, (S, S)>,
    params: &BranchPointOptimizationParams,
) -> Result<FindRootResult<S>, Report>
where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Send + Sync,
```

Trait bounds shrink from `ClockNode + ClockEdge` to `GraphNode + GraphEdge`. The statistics come from a pre-computed map rather than from node/edge payloads. This decouples the search from how statistics are stored on the graph.

`FindRootResult<S>` carries the winning `S` instead of `ClockSet`:

```rust
pub struct FindRootResult<S> {
    pub edge: Option<GraphEdgeKey>,
    pub split: f64,
    pub score: f64,
    pub stats: S,
}
```

## Statistics initialization

Each `RootStats` implementation needs a pre-search traversal to populate `edge_stats`:

- `ClockStats`: the existing clock regression pass (`clock_regression.rs`) already computes `to_parent`/`to_child` on each edge. Extract the edge statistics map from the regression result.
- `DivStats`: a new two-pass traversal (leaves-to-root, root-to-leaves) accumulating distance sums. Each edge gets `(to_parent: DivStats, to_child: DivStats)`.

The traversal for `DivStats` mirrors the clock regression's message-passing but only accumulates (count, d_sum, dsq_sum) from branch lengths. No dates, no clock model.

## Module placement

New top-level module `reroot/` in the `treetime` crate. Contains the generic algorithm extracted from `clock/find_best_root/`:

```
packages/treetime/src/reroot/
  mod.rs
  traits.rs            -- RootStats trait
  search.rs            -- find_best_root<S> (discrete node scan)
  split.rs             -- find_best_split<S> (1D refinement)
  cost_function.rs     -- EdgeCostFn<S> (argmin::CostFunction)
  method_brent.rs      -- moved from clock/find_best_root/
  method_golden_section.rs
  method_grid_search.rs
  params.rs            -- BranchPointOptimizationParams
  orchestrate.rs       -- reroot_in_place<S> (search + topology mutation)
  div_stats.rs         -- DivStats implementation
```

What stays in `clock/`:

```
packages/treetime/src/clock/
  reroot.rs            -- clock_reroot(): wraps reroot_in_place<ClockStats>,
                          select_root dispatch, clock edge fixup
  clock_regression.rs  -- estimate_clock_model_with_reroot_policy (unchanged API)
  find_best_root/      -- deleted (contents moved to reroot/)
```

What moves to `treetime-graph/`:

```
packages/treetime-graph/src/
  reroot.rs            -- unchanged (split_edge, apply_reroot_topology, remove_node_if_trivial)
  common_ancestor.rs   -- moved from clock/reroot.rs (LCA by path intersection, pure graph algorithm)
```

`common_ancestor` (`clock/reroot.rs:282-310`) uses only `GraphNode + GraphEdge` bounds with no domain logic. It belongs in the graph crate alongside other tree algorithms.

### Tip-based rerooting bypasses scoring

`--reroot-tips=A,B` finds MRCA via `common_ancestor`, roots at midpoint of parent edge (fixed split=0.5). No `RootStats` evaluation needed. This is a separate code path in `reroot/orchestrate.rs` that calls `common_ancestor` + `split_edge` + `apply_reroot_topology` directly, without `find_best_root<S>`.

### Downstream consumers after migration

| Consumer                     | Current call                                                    | After migration                                       |
| ---------------------------- | --------------------------------------------------------------- | ----------------------------------------------------- |
| `clock/pipeline.rs`          | `estimate_clock_model_with_reroot_policy`                       | unchanged (calls clock wrapper)                       |
| `timetree/pipeline.rs`       | `estimate_clock_model_with_reroot_policy` + `reroot_tree`       | unchanged                                             |
| `timetree/refinement.rs`     | `estimate_clock_model_with_reroot`                              | unchanged                                             |
| `commands/shared/reroot.rs`  | `RerootMethod`, `RerootSpec` from `clock/find_best_root/params` | from `clock/reroot` (dispatch enums stay clock-level) |
| `optimize/pipeline.rs` (new) | --                                                              | `reroot::orchestrate::reroot_in_place::<DivStats>`    |

Clock and timetree commands see no API change. Only `optimize` uses the generic path directly.

## Domain fixup after reroot

After `apply_reroot_topology` inverts edges, domain-specific data on those edges may need updating. Two cases:

- ClockStats: swap `to_parent`/`to_child` messages on inverted edges (existing `clock/reroot.rs:327-336`). Required because clock regression stores persistent per-edge statistics.
- DivStats: no-op. Statistics are ephemeral (computed for the search, discarded after). The partition fixup (`PartitionRerootOps::apply_reroot`) and post-reroot `update_marginal` handle the sequence reconstruction side.

The fixup is a callback provided by the caller of `reroot_in_place`, not a method on `RootStats`, since it operates on the graph's edge payload type (which `RootStats` is decoupled from):

```rust
pub fn reroot_in_place<S: RootStats>(
    graph: &mut Graph<N, E, D>,
    edge_stats: &BTreeMap<GraphEdgeKey, (S, S)>,
    params: &BranchPointOptimizationParams,
    reroot_params: &RerootParams,
    fixup: impl FnMut(&mut Graph<N, E, D>, &[GraphEdgeKey]) -> Result<(), Report>,
) -> Result<RerootResult<S>, Report>
```

Clock callers pass a closure that swaps ClockSet messages. Optimize callers pass a no-op.

### RerootParams factoring

The current `RerootParams` (`clock/reroot.rs:22-45`) bundles generic and clock-specific fields:

| Field                       | Generic or clock-specific | Disposition                                                                         |
| --------------------------- | ------------------------- | ----------------------------------------------------------------------------------- |
| `split_edge: bool`          | generic                   | stays in `RerootParams` (controls whether to split an edge or snap to nearest node) |
| `remove_trivial_root: bool` | generic                   | stays (topology cleanup after root move)                                            |
| `spec: RerootSpec`          | clock/command dispatch    | moves to clock wrapper (method -> objective mapping)                                |
| `objective: RootObjective`  | clock-specific            | moves to clock wrapper (EstimatedRate vs FixedRate scoring variant)                 |
| `force_positive_rate: bool` | clock-specific            | moves to clock wrapper (only meaningful with clock regression)                      |

The generic `RerootParams` retains only `split_edge` and `remove_trivial_root`. The clock wrapper (`clock/reroot.rs`) holds the clock-specific fields and maps them before calling `reroot_in_place`.

## Validation

### Behavioral equivalence

The refactoring must not change results for existing clock and timetree commands. Verify:

- Clock: run `clock --reroot=least-squares` on flu/h3n2/20, ebola, zika before and after. Root position (edge key + split fraction) and clock model parameters must match exactly.
- Timetree: run `timetree` on flu/h3n2/20 before and after. Timetree output (node times, branch lengths, clock model) must match within 1e-10.

### DivStats correctness

- Analytical: for a star tree (all tips at distance d from root), `DivStats.score()` = 0 at the true root (zero variance). For any other root position, score > 0.
- Synthetic: generate trees with known optimal root (midpoint), verify `find_best_root::<DivStats>` finds it.
- Cross-check: on datasets with dates, compare root position from `ClockStats` (MinDev mode, `force_positive_rate=false`) with `DivStats`. When dates carry no temporal signal, both should converge to the same root.

### Edge cases

- Single-tip tree: `DivStats.count = 1` at any root, score = 0 everywhere. Search returns current root (no improvement possible).
- Two-tip tree: optimal root at midpoint of the single edge. Score is a parabola in x with minimum at 0.5.
- All-zero branch lengths: `DivStats.score() = 0` everywhere. Search returns current root.

## Related

- [kb/proposals/optimize-pipeline-timetree-parity.md](optimize-pipeline-timetree-parity.md) -- parent proposal
- [kb/proposals/optimize-reroot-support.md](optimize-reroot-support.md) -- optimize-specific wiring using DivStats
- [kb/tickets/reroot-create-module-with-root-stats-trait.md](../tickets/reroot-create-module-with-root-stats-trait.md) -- create generic module and trait
- [kb/tickets/reroot-implement-div-stats-scoring.md](../tickets/reroot-implement-div-stats-scoring.md) -- DivStats implementation
- [kb/tickets/reroot-migrate-clock-to-generic-search.md](../tickets/reroot-migrate-clock-to-generic-search.md) -- migrate clock callers, delete old module
- [kb/features/optimize.md](../features/optimize.md) -- optimize feature inventory
