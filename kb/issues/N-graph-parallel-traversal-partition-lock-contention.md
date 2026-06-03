# Parallel traversal partition write locks serialize the frontier

Parallel BFS traversal callbacks in marginal and Fitch passes acquire `partition.write_arc()` for each partition on every node. Nodes in the same BFS frontier run concurrently via Rayon, but the write lock on the partition serializes them. The parallelism has no effect for partition-heavy workloads.

## Locations

- `packages/treetime/src/ancestral/marginal.rs` - `run_marginal_backward`, `run_marginal_forward` acquire `partition.write_arc()` per node
- `packages/treetime/src/ancestral/fitch.rs` - `run_fitch_forward` acquires `partition.nodes_mut()` per node

## Impact

Performance only. No correctness issue. BFS guarantees that frontier siblings are independent (all predecessors visited), so the write lock prevents data races. The contention degrades parallel traversal to serial for partition mutations.

## Possible approaches

- Per-node locking instead of per-partition (partition stores node data in a concurrent map)
- Batch partition updates after frontier completes (collect mutations, apply in bulk)
- Accept current behavior if profiling shows partition lock is not the bottleneck
