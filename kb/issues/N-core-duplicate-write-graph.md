# Duplicate `write_graph()` across commands

Three commands (`ancestral`, `optimize`, `prune`) each define a private `write_graph()` helper that serializes the final `GraphAncestral` to Newick and NEXUS. The implementations share the same generic signature and trait bounds, and differ only in the output filename prefix (`annotated_tree` for ancestral and optimize, `pruned_tree` for prune) and a commented-out JSON block. A shared helper parameterized by the filename prefix would eliminate the duplication.

## Locations

- Ancestral: [packages/treetime/src/commands/ancestral/run.rs#L208-L233](../../packages/treetime/src/commands/ancestral/run.rs#L208-L233) -- `annotated_tree.{nwk,nexus}`
- Optimize: [packages/treetime/src/commands/optimize/run.rs#L742-L767](../../packages/treetime/src/commands/optimize/run.rs#L742-L767) -- `annotated_tree.{nwk,nexus}`
- Prune: [packages/treetime/src/commands/prune/run.rs#L262-L281](../../packages/treetime/src/commands/prune/run.rs#L262-L281) -- `pruned_tree.{nwk,nexus}`

## Proposed placement

`packages/treetime/src/commands/shared/` or `packages/treetime-io/` (since this is a pure I/O helper). A single function like:

```rust
pub fn write_graph_files<N, E, D>(
  outdir: impl AsRef<Path>,
  stem: &str,
  graph: &Graph<N, E, D>,
) -> Result<(), Report>
```

with callers passing `"annotated_tree"` or `"pruned_tree"`.

## Related

- `packages/treetime/src/representation/algo/topology_cleanup/` -- shared module for cross-command graph operations. `write_graph()` does not belong there (it is pure I/O, not a graph algorithm) but the same consolidation pattern applies.
- `packages/treetime-io/` -- existing I/O crate; the shared helper naturally lives here.
