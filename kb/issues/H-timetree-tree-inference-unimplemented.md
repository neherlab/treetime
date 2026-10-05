# No tree inference from alignment

Every v1 command that reads a tree requires `--tree`: without it, the command stops with the required-arguments error of the command line. There is no code path that infers a tree from the alignment.

## Location

Each command declares `--tree` as `Option<PathBuf>`, but its argument conversion stops with the required-arguments error when the flag is missing:

- `timetree`: [`packages/app-commands/src/commands/timetree/args.rs#L157-L159`](../../packages/app-commands/src/commands/timetree/args.rs#L157-L159)
- `clock`: [`packages/app-commands/src/commands/clock/args.rs#L101-L108`](../../packages/app-commands/src/commands/clock/args.rs#L101-L108)
- `ancestral`: [`packages/app-commands/src/commands/ancestral/args.rs#L106`](../../packages/app-commands/src/commands/ancestral/args.rs#L106)
- `homoplasy`: [`packages/app-commands/src/commands/homoplasy/args.rs#L70-L72`](../../packages/app-commands/src/commands/homoplasy/args.rs#L70-L72)
- `mugration`: [`packages/app-commands/src/commands/mugration/args.rs#L84`](../../packages/app-commands/src/commands/mugration/args.rs#L84)
- `optimize`: [`packages/app-commands/src/commands/optimize/args.rs#L73`](../../packages/app-commands/src/commands/optimize/args.rs#L73)
- `prune`: [`packages/app-commands/src/commands/prune/args.rs#L45`](../../packages/app-commands/src/commands/prune/args.rs#L45)

## v0 behavior

v0 has no tree builder of its own. When no tree is provided, `fn tree_inference()` calls the external programs IQ-TREE, FastTree, and RAxML as subprocesses, in that order, and uses the first one that succeeds ([`packages/legacy/treetime/treetime/utils.py#L410-L452`](../../packages/legacy/treetime/treetime/utils.py#L410-L452)).

## Impact

Users who invoke `treetime timetree` without `--tree`, expecting tree inference from the alignment, get the required-arguments error. This blocks a standard v0 workflow.

## Fix

Implement tree inference from the alignment (a large feature tracked in [algo/unimplemented.md](../algo/unimplemented.md)), then make `--tree` optional again for the commands that can use it.

The design decisions are open: [kb/reports/feat-tree-infer.md](../reports/feat-tree-infer.md) surveys the methods and lists the options for the entry point (behavior without an input tree), initial topology, topology improvement, divergent data, branch support, update mode, and recombination. An internal builder diverges from v0, which only calls external programs, so each adopted option needs an approved entry in `kb/decisions/` before implementation.
