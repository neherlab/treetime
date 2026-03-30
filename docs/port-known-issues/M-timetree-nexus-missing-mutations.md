# Nexus output missing mutation annotations

v0 nexus annotates each branch with mutations:
`[&mutations="A55G,T93C",date=2003.84]`. v1 nexus only contains date:
`[&date="2000.68"]`.

## Status

Partially addressed by `PartitionBranchOps` trait promotion:

- `NodeAncestral.nwk_comments()` now emits mutations from the `mutations` field (populated by `annotate_branch_mutations()`)
- `NodeTimetree.nwk_comments()` now delegates to `NodeAncestral.nwk_comments()`, inheriting mutation support
- `annotate_branch_mutations()` called before Newick/Nexus output in ancestral and optimize commands

**Remaining**: timetree command does not yet call `annotate_branch_mutations()` because `PartitionBranchOps::edge_subs()` requires `&GraphAncestral` but the timetree operates on `Graph<NodeTimetree, EdgeTimetree, ()>`. The generic `edge_subs_from_graph()` method exists on `PartitionMarginalSparse` but cannot be used through the trait object interface. A timetree-specific annotation function is needed.

## Root cause

~~Three independent gaps:~~

1. ~~`NodeTimetree.nwk_comments()` builds the comment map from scratch, inserting only `date`.~~ **FIXED**: delegates to `NodeAncestral.nwk_comments()`

2. ~~`NodeAncestral.nwk_comments()` has a TODO stub.~~ **FIXED**: emits mutations from `NodeAncestral.mutations` field

3. `write_outputs()` in
   [`run.rs#L243`](../../packages/treetime/src/commands/timetree/run.rs#L243)
   receives `_partitions` (underscored, unused). The `annotate_branch_mutations()` function exists but requires `&GraphAncestral`, which the timetree command does not have.

## Related issues

- [Missing output files compared to v0](N-timetree-missing-output-files.md)
  lists `branch_mutations.txt` as another missing output
