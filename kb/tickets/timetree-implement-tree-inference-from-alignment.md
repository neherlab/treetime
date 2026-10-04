# Implement tree inference from alignment

The `timetree` and `clock` commands declare `--tree` as `Option<PathBuf>`, but their argument conversion stops with the required-arguments error when the flag is missing ([`packages/app-commands/src/commands/timetree/args.rs#L157-L159`](../../packages/app-commands/src/commands/timetree/args.rs#L157-L159), [`packages/app-commands/src/commands/clock/args.rs#L101-L108`](../../packages/app-commands/src/commands/clock/args.rs#L101-L108)). No code path infers a tree from the alignment.

v0 has no tree builder of its own. When no tree is provided, `fn tree_inference()` calls the external programs IQ-TREE, FastTree, and RAxML as subprocesses, in that order ([`packages/legacy/treetime/treetime/utils.py#L410-L452`](../../packages/legacy/treetime/treetime/utils.py#L410-L452)).

## Impact

Users who invoke `treetime timetree` without `--tree`, expecting tree inference from the alignment, get the required-arguments error. This blocks a standard v0 workflow.

## Fix

Implement tree inference from the alignment (a large feature tracked in [unimplemented.md](../algo/unimplemented.md)), then make `--tree` optional again for the commands that can use it.

This ticket is not ready for execution. The design decisions are open, and [kb/reports/feat-tree-infer.md](../reports/feat-tree-infer.md) lists them with options: behavior without an input tree, initial topology, topology improvement, divergent data, branch support, update mode, and recombination. An internal builder diverges from v0, so each adopted option needs an approved entry in `kb/decisions/`.

## Related issues

- Source: [kb/issues/H-timetree-tree-inference-unimplemented.md](../issues/H-timetree-tree-inference-unimplemented.md) -- delete after full resolution
- [unimplemented.md](../algo/unimplemented.md) -- tree inference algorithms listed as unimplemented
