# Dense marginal backward pass wipes indels from edge partitions

The `PartitionMarginalOps::process_node_backward()` implementation for `PartitionMarginalDense` at [packages/treetime/src/representation/partition/marginal_dense.rs#L211](../../packages/treetime/src/representation/partition/marginal_dense.rs#L211) creates a new `DenseEdgePartition::default()` for each edge and inserts it, replacing any existing edge data. This wipes the `indels: Vec<InDel>` field that was populated during Fitch reconstruction.

## Current behavior

Each call to `update_marginal()` for dense partitions:

1. Calls `process_node_backward()` on each node
2. For each parent edge, constructs `DenseEdgePartition::default()` (which has `indels: vec![]`)
3. Sets `msg_from_child` and `msg_to_parent` on the new edge data
4. Inserts the new edge data via `self.edges.insert(*edge_key, edge_data)`, overwriting any previous entry

As a result, indels set during Fitch reconstruction (`compress_sequences()`) are lost after the first `update_marginal()` call.

## Impact

Low in the current codebase. Indels for optimization are read from sparse partitions, which preserve indels across marginal updates. The `edge_indel_count()` method on `PartitionMarginalDense` always returns 0 after the first marginal pass, but the sparse partition provides the correct count.

If future code relies on dense partition indels for non-optimization purposes (e.g., ancestral sequence output with dense mode), the missing indels would be a correctness issue.

## Fix

Preserve the existing `indels` field when rebuilding edge partition data in the dense backward pass. The `msg_from_child` and `msg_to_parent` messages should be updated without creating a new `DenseEdgePartition` from scratch.
