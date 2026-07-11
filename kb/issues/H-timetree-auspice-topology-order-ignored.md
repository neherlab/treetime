# Timetree Auspice output ignores topology order

## Symptom

`--ladderize=ascending`, `--ladderize=descending`, and `--topology-order` have no effect on the Auspice v2 JSON output. The tree is written in the original traversal order regardless of the flag.

## Root cause

[`packages/treetime/src/commands/timetree/run.rs`](../../packages/treetime/src/commands/timetree/run.rs) builds a reordered graph and passes it to `write_tree_outputs`:

```rust
let plan = resolved.topology_order.plan(&output.graph)?;
let ordered = plan.ordered_graph(&output.graph)?;
// ...
write_tree_outputs(&ordered, &resolved.tree_outputs, &providers, Some(&auspice_ctx))?;
```

For Newick and Nexus output formats, [`packages/treetime-io/src/graph.rs`](../../packages/treetime-io/src/graph.rs) uses the `ordered` graph it receives. For `TreeWriteKind::Auspice`, it ignores the `graph` argument and delegates to the `AuspiceWriter` trait object:

```rust
TreeWriteKind::Auspice => {
  writer.write_auspice(path)?;
}
```

The `TimetreeAuspiceCtx` passed as `auspice_writer` holds a reference to the original `output.graph`. The Auspice JSON is therefore always written in the original node insertion order.

## Fix direction

Pass the graph supplied to `write_tree_outputs` into `AuspiceWriter::write_auspice`. Every tree format then serializes the same ordered topology, and the command-specific writer retains only the additional Auspice metadata it needs.
