# No tree inference from alignment

Every v1 command requires `--tree`: without it, the command stops with the required-arguments error of the command line. There is no code path that infers a tree from the alignment.

v0 infers a tree from the alignment using neighbor-joining or other methods when no tree is provided.

## Impact

Users who invoke `treetime timetree` without `--tree`, expecting tree inference from the alignment, get the required-arguments error. This blocks a standard v0 workflow.

## Fix

Implement tree inference from the alignment (a large feature tracked in [algo/unimplemented.md](../algo/unimplemented.md)), then make `--tree` optional again for the commands that can use it.

## Related tickets

- [kb/tickets/timetree-implement-tree-inference-from-alignment.md](../tickets/timetree-implement-tree-inference-from-alignment.md)
