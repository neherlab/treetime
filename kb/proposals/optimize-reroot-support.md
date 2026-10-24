# Reroot support for optimize command

Add `--reroot`, `--outgroup`, and `--keep-root` to the `optimize` command. Uses the date-free root search, which needs no dates and counts every tip ([kb/algo/reroot.md](../algo/reroot.md)), once on the branch lengths after the initial guess.

## Available methods

The optimize command has no date information. Only divergence-based rooting methods are valid:

| Method              | Valid                | Scoring                                                  |
| ------------------- | -------------------- | -------------------------------------------------------- |
| `MinDev`            | yes                  | rate-zero least squares -- minimize variance of root-to-tip distances |
| `Midpoint`          | yes                  | middle of the longest tip-to-tip path                    |
| `--outgroup=A,B` | yes                  | topology-only -- root at the middle of the branch above the MRCA of named tips |
| `LeastSquares`      | no -- requires dates | --                                                       |
| `Oldest`            | no -- requires dates | --                                                       |

Requesting a date-dependent method produces a CLI error.

## CLI args

Add to `TreetimeOptimizeArgs` (`commands/optimize/args.rs`):

```rust
#[cfg_attr(feature = "clap", clap(flatten))]
pub reroot: RerootArgs,

#[cfg_attr(feature = "clap", clap(long, conflicts_with_all = ["reroot", "outgroup"]))]
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
5. `reroot_date_free` (`reroot/date_free.rs`): the date-free reroot driver, on the branch lengths after the initial guess, with every tip counted
6. Apply reroot changes to partitions (`PartitionRerootOps::apply_reroot`)
7. Re-run `update_marginal` (root position changed, messages invalidated)
8. `run_optimize_loop` (main convergence loop)
9. Final `update_marginal`

The reroot uses the starting branch lengths, as the date-free methods do in every command ([kb/decisions/reroot-date-free-methods-root-once.md](../decisions/reroot-date-free-methods-root-once.md)), and the main loop then optimizes the branch lengths on the rooted tree.

When no reroot flags are given, steps 5-7 are skipped and behavior is identical to the current pipeline.

## Validation

- `optimize --reroot=min-dev` on flu/h3n2/20: root position changes, branch lengths re-optimized
- `optimize --outgroup=A,B` on flu/h3n2/20: root at MRCA of named tips
- `optimize --keep-root`: identical to current behavior
- Default (no flags): identical to current behavior
- `optimize --reroot=least-squares`: CLI error
- Root position from `optimize --reroot=min-dev` matches root from running `clock --reroot=min-dev` on the same tree

## Related

- [kb/proposals/optimize-pipeline-timetree-parity.md](optimize-pipeline-timetree-parity.md) -- parent proposal
- Status: implemented in `packages/app-commands/src/commands/optimize/` and `packages/treetime/src/optimize/pipeline.rs`
- Sibling pattern: `commands/timetree/args.rs`, `commands/clock/args.rs` (both have `RerootArgs`)
