# Shared output module mixes output planning with topology ordering

`commands/shared/output.rs` contains output kinds, command capability policy, CLI output fields, path derivation, and resolved output maps alongside topology-order arguments, target-order file loading, graph traversal, and private Newick payload types [`packages/treetime/src/commands/shared/output.rs`](../../packages/treetime/src/commands/shared/output.rs).

## Cohesion problem

Output-path changes and topology-order changes use different data and dependencies. Topology ordering pulls graph traversal and Newick parsing into a module whose primary responsibility is output selection and resolution. The source itself presents topology ordering as a separate concern.

The combined module obscures ownership: topology order influences tree materialization, but it is neither an output path nor a writer format. Keeping both in one file also makes the high-afferent output surface the import route for an independent graph operation.

## Required organization

Keep output descriptors, selection, defaults, and pure path planning in the output-planning module. Move topology-order parsing, target loading, and graph integration beside `treetime_graph::topology_order` or a focused application adapter. The output plan may reference a parsed ordering policy without owning how that policy is loaded or computed.

## Validation

- Output planning tests require no graph or Newick fixture.
- Topology-order tests cover target loading, duplicate labels, missing labels, and graph ordering independently.
- Every command produces the same ordered tree outputs after the move.
- Dependency analysis confirms output selection no longer imports graph traversal or Newick parsing solely for topology ordering.

## Related issues

- [M-command-output-ownership-is-scattered.md](M-command-output-ownership-is-scattered.md)
- [M-topology-order-duplicate-label-collapse.md](M-topology-order-duplicate-label-collapse.md)
