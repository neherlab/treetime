# Sparse reconstruction skips non-root nodes without exactly one parent

## Summary

Sparse marginal reconstruction treats a failed parent lookup as "do not emit this node". A non-root node with zero or several parents, which a graph with reticulations can contain, is silently left without a reconstructed sequence instead of producing an error.

## Current state

`PartitionMarginalSparse::advance_node_state()` [packages/treetime/src/partition/marginal/sparse/partition.rs#L129](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L129) reads the parent of every non-root node with `get_exactly_one(&node.parent_keys).ok()?`. The function returns `Option<()>`, where `None` already means "the node is not emitted" (a leaf when leaves are excluded). The `.ok()?` conversion merges a second meaning into the same `None`: the graph contradicts the single-parent assumption of the forward message.

Callers in `ancestral_reconstruction` and the timetree reconstruction therefore cannot distinguish an intentionally skipped leaf from a node whose parent structure the sparse forward pass does not support.

## Required contract

- Node emission and graph-shape failures travel on separate channels: the reconstruction step returns `Result<Option<Seq>, Report>`.
- A non-root node without exactly one parent produces an error that names the node and its parent count, until the sparse forward pass defines how several parent messages combine.
- Tree inputs keep their current output.

## Related issues

- [H-graph-capability-contracts-silently-discard-state.md](H-graph-capability-contracts-silently-discard-state.md)
