# Partition read views are pull accessors, not value maps

Inference stages return per-node and per-edge values, but the consumers of a completed reconstruction do not read those values directly. They go through two accessor traits and a dynamic view: [`packages/treetime/src/partition/traits.rs#L88`](../../packages/treetime/src/partition/traits.rs#L88) declares `trait PartitionBranchOps` (sequence length, edge substitutions, edge indels, root sequence, node sequence, effective length) and [`packages/treetime/src/partition/traits.rs#L131`](../../packages/treetime/src/partition/traits.rs#L131) declares `trait PartitionOptimizeOps` (per-edge likelihood contribution, indel count). The optimizer builds `Vec<&dyn PartitionOptimizeOps>` at [`packages/treetime/src/optimize/run_loop.rs#L78`](../../packages/treetime/src/optimize/run_loop.rs#L78) and every consumer calls back per edge.

Each accessor returns a backend-neutral value (`Vec<Sub>`, `Seq`, `usize`, `OptimizationContribution`), so nothing in the results requires dynamic dispatch: one precomputed map per quantity would serve the same consumers. The pull design keeps the whole reconstruction alive throughout every consumer, hides which quantities a consumer actually reads, and makes each read a virtual call inside per-edge loops.

## Required contract

- A stage that produces per-edge or per-node quantities publishes them as explicit maps keyed by `GraphEdgeKey` or `GraphNodeKey`.
- Output writers and the branch-length optimizer read those maps and hold no reference to a partition or reconstruction.
- Dynamic dispatch over representations remains only where a genuine backend difference exists, not for value reads.

## Dependency

The inner optimizer functions take `&[&dyn PartitionOptimizeOps]`, so the maps must exist before the traits can be deleted. The shape of the map set also depends on whether a single run can hold both dense and sparse partitions at once (see [N-partition-mixed-dense-sparse-view-unreachable.md](N-partition-mixed-dense-sparse-view-unreachable.md)): a homogeneous run needs one map per quantity, a mixed run needs a defined combination rule per quantity.

## Validation

- Optimize, ancestral, timetree, and prune outputs stay byte-identical across the full smoke matrix.
- No production type implements a partition read trait after the migration.
- Per-edge optimizer loops read from a map, not through dynamic dispatch.
