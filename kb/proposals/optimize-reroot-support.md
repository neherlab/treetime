# Reroot support for optimize command

Add `--reroot`, `--reroot-tips`, and `--keep-root` to the `optimize` command. Uses the generic root search architecture with `DivStats` scoring (divergence-only, no dates required). Implements a two-phase branch-length optimization pattern around the reroot step, matching timetree's pipeline structure.

Depends on [reroot-generic-scoring-architecture.md](reroot-generic-scoring-architecture.md).

## Available methods

The optimize command has no date information. Only divergence-based rooting methods are valid:

| Method              | Valid                | Scoring                                                  |
| ------------------- | -------------------- | -------------------------------------------------------- |
| `MinDev`            | yes                  | `DivStats` -- minimize variance of root-to-tip distances |
| `--reroot-tips=A,B` | yes                  | topology-only -- root at MRCA of named tips              |
| `LeastSquares`      | no -- requires dates | --                                                       |
| `Oldest`            | no -- requires dates | --                                                       |
| `ClockFilter`       | no -- requires dates | --                                                       |

Requesting a date-dependent method produces a CLI error.

## CLI args

Add to `TreetimeOptimizeArgs` (`commands/optimize/args.rs`):

```rust
#[cfg_attr(feature = "clap", clap(flatten))]
pub reroot: RerootArgs,

#[cfg_attr(feature = "clap", clap(long, conflicts_with_all = ["reroot", "reroot_tips"]))]
pub keep_root: bool,
```

`RerootArgs` is the existing shared struct from `commands/shared/reroot.rs` (already used by clock and timetree).

## Pipeline changes

Current pipeline (`optimize/pipeline.rs`):

1. Create partition
2. Initialize + update marginal
3. Normalize rates
4. Apply initial guess
5. `run_optimize_loop`
6. Final `update_marginal`

New pipeline when `--reroot` is set and `--keep-root` is not:

1. Create partition
2. Initialize + update marginal
3. Normalize rates
4. Apply initial guess
5. Pre-reroot BL optimization pass (single `run_optimize_mixed` with damping 0.75)
6. Compute `DivStats` edge statistics (two-pass traversal over branch lengths)
7. `reroot_in_place::<DivStats>` with no-op fixup
8. Apply reroot changes to partitions (`PartitionRerootOps::apply_reroot`)
9. Re-run `update_marginal` (root position changed, messages invalidated)
10. `run_optimize_loop` (main convergence loop)
11. Final `update_marginal`

Steps 5-9 mirror timetree's two-phase pattern (`timetree/pipeline.rs:172-209`): optimize branch lengths first (so the root search uses reasonable distances), then reroot, then re-optimize (because the root position changes the optimization landscape).

When no reroot flags are given, steps 5-9 are skipped and behavior is identical to the current pipeline.

## Validation

- `optimize --reroot=min-dev` on flu/h3n2/20: root position changes, branch lengths re-optimized
- `optimize --reroot-tips=A,B` on flu/h3n2/20: root at MRCA of named tips
- `optimize --keep-root`: identical to current behavior
- Default (no flags): identical to current behavior
- `optimize --reroot=least-squares`: CLI error
- Root position from `optimize --reroot=min-dev` matches root from running `clock --reroot=min-dev` on the same tree

## Related

- [kb/proposals/reroot-generic-scoring-architecture.md](reroot-generic-scoring-architecture.md) -- prerequisite
- [kb/proposals/optimize-pipeline-timetree-parity.md](optimize-pipeline-timetree-parity.md) -- parent proposal
- Status: implemented in `packages/treetime/src/commands/optimize/` and `packages/treetime/src/optimize/pipeline.rs`
- Sibling pattern: `commands/timetree/args.rs`, `commands/clock/args.rs` (both have `RerootArgs`)
