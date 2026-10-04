# No tree inference from alignment

Every v1 command requires `--tree`: without it, the command stops with the required-arguments error of the command line. There is no code path that infers a tree from the alignment.

v0 has no tree builder of its own. When no tree is provided, `fn tree_inference()` calls the external programs IQ-TREE, FastTree, and RAxML as subprocesses, in that order, and uses the first one that succeeds ([`packages/legacy/treetime/treetime/utils.py#L410-L452`](../../packages/legacy/treetime/treetime/utils.py#L410-L452)).

## Impact

Users who invoke `treetime timetree` without `--tree`, expecting tree inference from the alignment, get the required-arguments error. This blocks a standard v0 workflow.

## Fix

Implement tree inference from the alignment (a large feature tracked in [algo/unimplemented.md](../algo/unimplemented.md)), then make `--tree` optional again for the commands that can use it. The design decisions are open: [kb/reports/feat-tree-infer.md](../reports/feat-tree-infer.md) surveys the methods and lists the options for the entry point, initial topology, topology improvement, divergent data, branch support, update mode, and recombination. An internal builder diverges from v0, which only calls external programs.

## Related tickets

- [kb/tickets/timetree-implement-tree-inference-from-alignment.md](../tickets/timetree-implement-tree-inference-from-alignment.md)
