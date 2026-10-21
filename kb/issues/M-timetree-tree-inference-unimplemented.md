# No tree inference from alignment

Every v1 command that reads a tree requires `--tree`: without it, the command stops with the required-arguments error of the command line. There is no code path that infers a tree from the alignment.

## Location

Each command declares `--tree` as `Option<PathBuf>`, but its argument conversion stops with the required-arguments error when the flag is missing:

- `timetree`: [`packages/app-commands/src/commands/timetree/args.rs#L64-L67`](../../packages/app-commands/src/commands/timetree/args.rs#L64-L67)
- `clock`: [`packages/app-commands/src/commands/clock/args.rs#L42-L49`](../../packages/app-commands/src/commands/clock/args.rs#L42-L49)
- `ancestral`: [`packages/app-commands/src/commands/ancestral/args.rs#L64-L67`](../../packages/app-commands/src/commands/ancestral/args.rs#L64-L67)
- `homoplasy`: [`packages/app-commands/src/commands/homoplasy/args.rs#L55-L58`](../../packages/app-commands/src/commands/homoplasy/args.rs#L55-L58)
- `mugration`: [`packages/app-commands/src/commands/mugration/args.rs#L41-L60`](../../packages/app-commands/src/commands/mugration/args.rs#L41-L60)
- `optimize`: [`packages/app-commands/src/commands/optimize/args.rs#L74-L77`](../../packages/app-commands/src/commands/optimize/args.rs#L74-L77)
- `prune`: [`packages/app-commands/src/commands/prune/args.rs#L36-L39`](../../packages/app-commands/src/commands/prune/args.rs#L36-L39)

## v0 behavior

v0 has no tree builder of its own. When no tree is provided, `fn tree_inference()` calls the external programs IQ-TREE, FastTree, and RAxML as subprocesses, in that order, and uses the first one that succeeds ([`packages/legacy/treetime/treetime/utils.py#L410-L452`](../../packages/legacy/treetime/treetime/utils.py#L410-L452)).

## Impact

Users who invoke `treetime timetree` without `--tree`, expecting tree inference from the alignment, get the required-arguments error. This blocks a standard v0 workflow.

## Fix

Implement tree inference from the alignment (a large feature tracked in [kb/algo/unimplemented.md](../algo/unimplemented.md)), then accept a missing `--tree` in the commands that can use an inferred tree.

[kb/reports/feat-tree-infer.md](../reports/feat-tree-infer.md) surveys the methods and lists the open design decisions: the entry point, initial topology, topology improvement, divergent data, branch support, update mode, and recombination. An internal builder diverges from v0, which only calls external programs, so each adopted option needs an approved entry in `kb/decisions/` before implementation.

> [!IMPORTANT]
> **Decision required.** v1 has no answer for a run without `--tree`, and the v0 behavior depends on external programs. Options for the entry point:
>
> - **Keep the error**: the user builds the tree with an external program first, for example `augur tree`
> - **Call an external program** as v0 does: run IQ-TREE, FastTree, or RAxML as a subprocess when one is installed. The report finds GPL or AGPL licenses for all production ML builders, so TreeTime can run them only as separate processes
> - **Build the tree internally**: a native builder for low-divergence data, with an error that names an external builder for divergent data
>
> The report lists six further decisions for an internal builder. Evidence: the v0 code above and the survey in [kb/reports/feat-tree-infer.md](../reports/feat-tree-infer.md)
