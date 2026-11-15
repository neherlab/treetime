# optimize cannot fit branch lengths on graphs of segment trees

`optimize` takes one tree. On a graph of segment trees ([kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md)) it is the first command that combines the segments: an edge shared by several segments has one length, and its objective is the sum of the contributions of every partition that the edge carries. Ancestral reconstruction treats segments independently once the branch lengths are fixed; branch-length optimization does not.

## Work needed

- **Partitions per segment**: `PartitionSet` built with the `Inheritance` of each segment, as `ancestral` does ([packages/treetime/src/ancestral/plan.rs](../../packages/treetime/src/ancestral/plan.rs))
- **Edge objective**: the contributions of an edge summed over the partitions that carry it; an edge that no partition carries keeps its length ([packages/treetime/src/partition/optimize/contribution.rs](../../packages/treetime/src/partition/optimize/contribution.rs))
- **Initial guess and rate normalization** over the carried edges and the partitions of all segments
- **Topology edits and rerooting** (polytomy merging, reversion hoisting, zero-length collapse, reroot) are off on a segment graph until the edge scopes have rules for them ([kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md#topology-operations))
- **Command interface**: `--segment-trees` and `--segment-alignments` as in `ancestral`
- **Pass-through nodes**: without them, a segment that passes through a node of another segment has its own edge beside the shared path, so the lengths of the two are fitted separately ([kb/reports/2026-11-15-treeknit-output-to-named-segment-trees.md](../reports/2026-11-15-treeknit-output-to-named-segment-trees.md))

## Validation

- Every segment carries every edge: the result equals the multi-partition `optimize` of one tree
- No edge is shared: each segment's lengths equal a separate `optimize` of its tree
- Toy reassortment graph: a shared edge's length maximizes the summed likelihood of its segments
