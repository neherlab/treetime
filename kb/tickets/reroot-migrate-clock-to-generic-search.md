# Migrate clock rerooting to generic search via ClockStats

Make the existing clock rerooting use the generic `reroot/` search infrastructure. Wrap `ClockSet` with a `RootStats` implementation (`ClockStats`). Update `clock/reroot.rs` to call `reroot::find_best_root::<ClockStats>` instead of the old `clock/find_best_root::find_best_root`. Delete `clock/find_best_root/` after migration.

Builds on the generic `reroot/` module (`RootStats` trait, `EdgeCostFn<S>`, generic search) already implemented in `packages/treetime/src/reroot/`.

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

### Set the reroot split-position optimizer default to Brent (v0 parity)

The timetree and clock reroot paths currently default the per-edge split-position optimizer to an 11-point grid. `timetree/pipeline.rs:134` builds `BranchPointOptimizationParams::default()`, whose `#[default]` variant is `Grid(GridSearchParams { n_points: 11 })` (`clock/find_best_root/params.rs`), and threads it through `timetree/refinement.rs` into `timetree/optimization/reroot.rs::reroot_tree` and `estimate_clock_model_with_reroot_policy`. The same `::default()` is used at two sites in `packages/treetime/examples/timetree_validation.rs`.

v0 refines the root position with Brent: `TreeRegression._optimal_root_along_branch` (`packages/legacy/treetime/treetime/treeregression.py`) evaluates a 6-point grid only to bracket, then calls `scipy.optimize.minimize_scalar(method='bounded', xatol=1e-6)`. An 11-point grid resolves the root to only ~0.1 branch fractions, diverging from v0. The optimize command's reroot already uses Brent (`optimize/pipeline.rs` builds `reroot::params::BranchPointOptimizationParams::brent()`), so the two reroot paths currently use different defaults; this migration is the point to make timetree and clock consistent.

Change: once the reroot path is on the generic module, construct the optimizer as `reroot::params::BranchPointOptimizationParams::brent()` instead of `::default()` at `timetree/pipeline.rs:134` and the two `timetree_validation.rs` sites (or wherever the default is centralized after migration).

This is a v0-parity behavior change for the clock and timetree commands. It intentionally shifts the rerooted split position relative to pre-migration v1 (finer resolution), so it cannot coexist with the "identical before and after" equivalence checks below as written -- see the note under Tests. Obtain team sign-off before implementing, since it changes numerical output.

## Tests

### Behavioral equivalence (zero regressions)

These checks isolate the migration mechanics from the Brent-default change above. Construct the split optimizer explicitly as `Grid(GridSearchParams { n_points: 11 })` on both sides so the root position is comparable; the default flip to Brent is verified separately (see below). Migrating with the new Brent default at the same time would shift the root and make these "identical before and after" comparisons fail for the wrong reason.

- `clock --reroot=least-squares` on flu/h3n2/20: root position (edge key + split fraction) and clock model parameters identical before and after
- `clock --reroot=min-dev` on flu/h3n2/20: identical
- `clock --reroot-tips=A,B` on flu/h3n2/20: identical
- `timetree` on flu/h3n2/20: timetree output identical before and after
- All existing tests in `clock/__tests__/test_reroot.rs` pass unchanged
- All existing tests in `timetree/optimization/__tests__/test_reroot.rs` pass unchanged

### Reroot default method (v0 parity)

- With the Brent default, `clock --reroot=min-dev` and `timetree` root split positions on flu/h3n2/20 match v0 `treetime` (which uses `minimize_scalar` Brent refinement) within numerical tolerance, where the old 11-point grid did not

### ClockStats unit tests

- `ClockSet` `RootStats::score()` matches `ClockSet::chisq()` for several known ClockSet values
- `ClockSet` `RootStats::propagate()` matches `ClockSet::propagate_averages()` for several inputs

## Related issues

- Proposal: [kb/proposals/reroot-generic-scoring-architecture.md](../proposals/reroot-generic-scoring-architecture.md)
- Builds on: `packages/treetime/src/reroot/` (generic module with `RootStats` and `DivStats`, already implemented)
- Fix during migration: [kb/issues/M-clock-mindev-wrong-objective.md](../issues/M-clock-mindev-wrong-objective.md) -- change MinDev from `EstimatedRate` to `FixedRate(0.0)`
