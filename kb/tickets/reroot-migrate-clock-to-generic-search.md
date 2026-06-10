# Migrate clock rerooting to generic search via ClockStats

Make the existing clock rerooting use the generic `reroot/` search infrastructure. Wrap `ClockSet` with a `RootStats` implementation (`ClockStats`). Update `clock/reroot.rs` to call `reroot::find_best_root::<ClockStats>` instead of the old `clock/find_best_root::find_best_root`. Delete `clock/find_best_root/` after migration.

Depends on [reroot-create-module-with-root-stats-trait.md](reroot-create-module-with-root-stats-trait.md).

## Scope

### ClockStats (`clock/clock_stats.rs` or impl on `ClockSet`)

Implement `RootStats` for `ClockSet` (or a newtype `ClockStats`):

- `leaf` -> `ClockSet::leaf_contribution_to_parent`
- `propagate` -> `ClockSet::propagate_averages`
- `score` -> `ClockSet::chisq()`

Decision: newtype vs direct impl. Direct impl on `ClockSet` is simpler (no wrapping/unwrapping). Newtype is cleaner if `ClockSet` should remain independent of the `reroot` crate. Recommend direct impl since `ClockSet` and `reroot` are in the same crate (`treetime`).

### Update clock/reroot.rs

- `select_root`: replace calls to `clock::find_best_root::find_best_root` with `reroot::search::find_best_root::<ClockSet>`
- `reroot_in_place`: replace internal search with `reroot::orchestrate::reroot_in_place::<ClockSet>`, passing a fixup closure that swaps `to_parent`/`to_child` on inverted edges
- Extract `ClockSet` edge statistics from the existing clock regression data into the `BTreeMap<GraphEdgeKey, (ClockSet, ClockSet)>` format expected by the generic search
- Delete `BranchPointCostFunction` -- `EdgeCostFn<ClockSet>` replaces it

### Update downstream imports

All files that import from `clock::find_best_root::params` need updating:

- `clock/clock_regression.rs` -- `BranchPointOptimizationParams` from `reroot::params`
- `clock/pipeline.rs` -- same
- `commands/clock/args.rs` -- `BrentParams`, `GoldenSectionParams`, `GridSearchParams`, `OptimizationMethod` from `reroot::params`
- `commands/clock/run.rs` -- `BranchPointOptimizationParams`, `OptimizationMethod` from `reroot::params`
- `commands/shared/reroot.rs` -- `RerootMethod`, `RerootSpec` stay in `clock/` (not moved)
- `timetree/optimization/reroot.rs` -- `BranchPointOptimizationParams`, `RerootSpec`
- `timetree/pipeline.rs` -- same
- `timetree/refinement.rs` -- `BranchPointOptimizationParams`

### Delete clock/find_best_root/

After all callers are migrated, delete the entire `clock/find_best_root/` directory:

- `find_best_root.rs` -- replaced by `reroot/search.rs`
- `find_best_split.rs` -- replaced by `reroot/split.rs`
- `cost_function.rs` -- replaced by `reroot/cost_function.rs`
- `method_brent.rs` -- moved to `reroot/`
- `method_golden_section.rs` -- moved to `reroot/`
- `method_grid_search.rs` -- moved to `reroot/`
- `params.rs` -- split between `reroot/params.rs` and `clock/reroot.rs`
- `mod.rs` -- deleted

## Tests

### Behavioral equivalence (zero regressions)

- `clock --reroot=least-squares` on flu/h3n2/20: root position (edge key + split fraction) and clock model parameters identical before and after
- `clock --reroot=min-dev` on flu/h3n2/20: identical
- `clock --reroot-tips=A,B` on flu/h3n2/20: identical
- `timetree` on flu/h3n2/20: timetree output identical before and after
- All existing tests in `clock/__tests__/test_reroot.rs` pass unchanged
- All existing tests in `timetree/optimization/__tests__/test_reroot.rs` pass unchanged

### ClockStats unit tests

- `ClockSet` `RootStats::score()` matches `ClockSet::chisq()` for several known ClockSet values
- `ClockSet` `RootStats::propagate()` matches `ClockSet::propagate_averages()` for several inputs

## Related issues

- Proposal: [kb/proposals/reroot-generic-scoring-architecture.md](../proposals/reroot-generic-scoring-architecture.md)
- Depends on: [kb/tickets/reroot-create-module-with-root-stats-trait.md](reroot-create-module-with-root-stats-trait.md)
- Independent of: [kb/tickets/reroot-implement-div-stats-scoring.md](reroot-implement-div-stats-scoring.md)
- Fix during migration: [kb/issues/M-clock-mindev-wrong-objective.md](../issues/M-clock-mindev-wrong-objective.md) -- change MinDev from `EstimatedRate` to `FixedRate(0.0)`
