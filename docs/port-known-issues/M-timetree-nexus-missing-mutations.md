# Nexus output missing mutation annotations

v0 nexus annotates each branch with mutations:
`[&mutations="A55G,T93C",date=2003.84]`. v1 nexus only contains date:
`[&date="2000.68"]`.

## Root cause

Three independent gaps:

1. `NodeTimetree.nwk_comments()` at
   [`timetree.rs#L129`](../../packages/treetime/src/representation/payload/timetree.rs#L129)
   builds the comment map from scratch, inserting only `date`. It does not
   delegate to `NodeAncestral.nwk_comments()`.

2. `NodeAncestral.nwk_comments()` at
   [`ancestral.rs#L35`](../../packages/treetime/src/representation/payload/ancestral.rs#L35)
   has a TODO stub: `let mutations: String = "".to_owned(); // TODO: fill mutations`.

3. `write_outputs()` in
   [`run.rs#L243`](../../packages/treetime/src/commands/timetree/run.rs#L243)
   receives `_partitions` (underscored, unused). Branch substitutions live in
   `SparseEdgePartition.subs` (partition layer), structurally inaccessible to
   the nexus writer.

## Related issues

- [Missing output files compared to v0](N-timetree-missing-output-files.md)
  lists `branch_mutations.txt` as another missing output
