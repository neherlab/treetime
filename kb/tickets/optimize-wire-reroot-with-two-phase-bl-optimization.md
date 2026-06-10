# Wire reroot into optimize with two-phase BL optimization

Add `--reroot`, `--reroot-tips`, and `--keep-root` to the optimize command. Implement a two-phase branch-length optimization pattern around the reroot step. Uses the generic root search architecture with `DivStats` scoring.

Depends on [reroot-create-module-with-root-stats-trait.md](reroot-create-module-with-root-stats-trait.md) and [reroot-implement-div-stats-scoring.md](reroot-implement-div-stats-scoring.md).

## Scope

### CLI args (`commands/optimize/args.rs`)

- Add `pub reroot: RerootArgs` (flatten)
- Add `pub keep_root: bool` with `conflicts_with_all = ["reroot", "reroot_tips"]`
- Runtime validation: reject `LeastSquares`, `Oldest`, `ClockFilter` (require dates)

### Pipeline (`optimize/pipeline.rs`)

When `--reroot` is set and `--keep-root` is not, insert steps between initial guess and the main loop:

1. Pre-reroot BL optimization: single `run_optimize_mixed` pass with damping 0.75
2. Compute `DivStats` edge statistics via two-pass traversal
3. `reroot::orchestrate::reroot_in_place::<DivStats>` with no-op fixup
4. `PartitionRerootOps::apply_reroot` on each partition
5. `update_marginal` to refresh messages after root move

When no reroot flags are given, these steps are skipped.

### Method dispatch

`RerootSpec::Method(MinDev)` -> call `reroot_in_place::<DivStats>` with DivStats edge statistics
`RerootSpec::Tips(names)` -> find MRCA via `common_ancestor`, root at midpoint of parent edge, no scoring needed
Date-dependent methods -> error before pipeline starts

## Locations

- CLI args: `packages/treetime/src/commands/optimize/args.rs`
- Pipeline: `packages/treetime/src/optimize/pipeline.rs`
- Command runner: `packages/treetime/src/commands/optimize/run.rs`
- Shared reroot args: `packages/treetime/src/commands/shared/reroot.rs`
- Generic search (from prerequisite): `packages/treetime/src/reroot/`
- Timetree reference pattern: `packages/treetime/src/timetree/pipeline.rs:172-209`

## Tests

- `optimize --reroot=min-dev` on flu/h3n2/20: root changes, BLs re-optimized
- `optimize --reroot-tips=TIP_A,TIP_B` on flu/h3n2/20: root at MRCA
- `optimize --keep-root`: identical to current behavior
- Default (no flags): identical to current behavior
- `optimize --reroot=least-squares`: error
- `optimize --reroot-tips=NONEXISTENT`: error
- Root position from `optimize --reroot=min-dev` compared with `clock --reroot=min-dev` on same input

## Related issues

- Proposal: [kb/proposals/optimize-reroot-support.md](../proposals/optimize-reroot-support.md)
- Depends on: [kb/tickets/reroot-create-module-with-root-stats-trait.md](reroot-create-module-with-root-stats-trait.md)
- Depends on: [kb/tickets/reroot-implement-div-stats-scoring.md](reroot-implement-div-stats-scoring.md)
- Parent: [kb/proposals/optimize-pipeline-timetree-parity.md](../proposals/optimize-pipeline-timetree-parity.md)
- Sibling pattern: `commands/timetree/args.rs`, `commands/clock/args.rs` (both have `RerootArgs`)
