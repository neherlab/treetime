# Create reroot module with RootStats trait and generic search

Extract the root search algorithm from `clock/find_best_root/` into a new top-level `reroot/` module. Define the `RootStats` trait that makes scoring pluggable. Implement generic `find_best_root<S>` and `find_best_split<S>`. Move the 1D optimization methods.

This ticket creates the generic infrastructure. Concrete `RootStats` implementations (`DivStats`, `ClockStats`) and caller migration are separate tickets.

## Scope

### RootStats trait (`reroot/traits.rs`)

```rust
pub trait RootStats: Clone + Default + Add<Output = Self> + Send + Sync {
    fn leaf(time: Option<f64>, branch_length: f64, variance: f64) -> Self;
    fn propagate(&self, branch_length: f64, variance: f64) -> Self;
    fn score(&self) -> f64;
}
```

### EdgeCostFn<S> (`reroot/cost_function.rs`)

Generic replacement for `BranchPointCostFunction`. Holds `to_parent: S`, `to_child: S`, branch length, variance, leaf flag. Implements `argmin::CostFunction<Param = f64, Output = f64>`.

### Generic search (`reroot/search.rs`, `reroot/split.rs`)

- `find_best_root<S: RootStats>`: iterate nodes, evaluate `S::score()` at each, return best. Trait bounds: `N: GraphNode + Named, E: GraphEdge` (no `ClockNode`/`ClockEdge`). Statistics come from a `BTreeMap<GraphEdgeKey, (S, S)>` parameter, not from node/edge payloads.
- `find_best_split<S: RootStats>`: optimize split position on candidate edge via 1D minimizer. Dispatches to method_brent/golden_section/grid_search.
- `FindRootResult<S>`: carries winning `S` stats, edge key, split fraction, score.

### Orchestration (`reroot/orchestrate.rs`)

`reroot_in_place<S: RootStats>`: search for best root, apply topology mutation (`treetime-graph` reroot), call fixup callback for domain-specific edge data updates. Takes `fixup: impl FnMut(...)` parameter.

### Move method files

Move from `clock/find_best_root/` to `reroot/`:

- `method_brent.rs` -- update imports, parameterize over `EdgeCostFn<S>` instead of `BranchPointCostFunction`
- `method_golden_section.rs` -- same
- `method_grid_search.rs` -- same

These already use `argmin::CostFunction` trait -- only the concrete type changes.

### Move params

`BranchPointOptimizationParams`, `GridSearchParams`, `BrentParams`, `GoldenSectionParams`, `OptimizationMethod` move to `reroot/params.rs`.

`RerootMethod`, `RerootSpec`, `RootObjective` stay in `clock/` -- they are dispatch/policy types, not search algorithm types.

## Locations

New files:

- `packages/treetime/src/reroot/mod.rs`
- `packages/treetime/src/reroot/traits.rs`
- `packages/treetime/src/reroot/search.rs`
- `packages/treetime/src/reroot/split.rs`
- `packages/treetime/src/reroot/cost_function.rs`
- `packages/treetime/src/reroot/orchestrate.rs`
- `packages/treetime/src/reroot/params.rs`
- `packages/treetime/src/reroot/method_brent.rs`
- `packages/treetime/src/reroot/method_golden_section.rs`
- `packages/treetime/src/reroot/method_grid_search.rs`

Files to keep temporarily (callers not yet migrated):

- `packages/treetime/src/clock/find_best_root/` -- still used by clock callers until migration ticket

## Tests

- `RootStats` trait: compile-time verification with a trivial test implementation
- `EdgeCostFn<S>`: unit test with a mock `RootStats` that returns known scores
- `find_best_root<S>`: unit test on a small synthetic tree with a trivial scoring function
- Method files: existing tests continue to pass (only imports change)

## Related issues

- Proposal: [kb/proposals/reroot-generic-scoring-architecture.md](../proposals/reroot-generic-scoring-architecture.md)
- Blocks: [kb/tickets/reroot-implement-div-stats-scoring.md](reroot-implement-div-stats-scoring.md)
- Blocks: [kb/tickets/reroot-migrate-clock-to-generic-search.md](reroot-migrate-clock-to-generic-search.md)
