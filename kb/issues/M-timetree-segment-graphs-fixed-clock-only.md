# Timetree on graphs of segment trees runs only at a fixed clock rate

`timetree --segment-trees --segment-alignments --clock-rate` dates a graph of segment trees merged by node name ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)). It optimizes the branch lengths jointly over the segments, then runs the time passes with several parents per node as loopy message passing, with refinement rounds that feed the clock branch lengths back into the marginal reconstruction ([kb/algo/timetree.md](../algo/timetree.md#graphs-with-several-parents)). Validation: [kb/reports/2026-11-15-segment-validation-timetree.md](../reports/2026-11-15-segment-validation-timetree.md).

`fn run_segment_timetree` [packages/treetime/src/timetree/segments.rs](../../packages/treetime/src/timetree/segments.rs) is a reduced pipeline. Compared to the tree pipeline it lacks:

- **Clock rate estimation**: `--clock-rate` is required, and `--clock-std-dev` is rejected. The root-to-tip regression needs one divergence per tip, which a tip with several root paths does not have; the marginal likelihood probe of the rate should carry over
- **Coalescent prior**: rejected. The expected lineage count of a segment graph is per segment, and the prior's merger events are the branching events of each segment
- **Rooting**: the roots of the graph are kept; `--reroot` and `--outgroup` are rejected. Rerooting inverts edges, which has no meaning at a node with several parents
- **Polytomy resolution, short-branch collapse, clock filter, tip pruning, relaxed clock**: not run, or rejected
- **Outputs**: reconstructed sequences, augur node data, the clock model and the root-to-tip table are not written; the node dates table covers every node of the graph

Known limitations of the passes:

- **Hidden nodes**: a node with one child in every segment (a reassortment node given as a pass-through node of every segment, or a merger of lineages with disjoint segments) is placed at its child with a narrow region, because optimize holds the edge below it at length 0, and it pulls the neighbouring deep nodes. Tree inference never produces such nodes, so real inputs lack them; an input that has them could be rejected, or the nodes contracted
- **Zero-length edges** (decided: kept): the flu graph has 391 nodes with a zero-length parent edge among 801 with several parents. TreeKnit expects collapsed input trees and resolves their polytomies with splits from compatible other segments, so such an edge is a branch that another segment supports, not an arbitrary resolution. The graph keeps them and does not contract them; inputs to TreeKnit should be collapsed trees
- **Cross-parent messages**: off by default (decided: little information in practice, slower convergence); `--cross-parent-messages` switches them on. With them, the sweep limit of 20 is reached on the flu graph at a relative change of 2 x 10^-4
- **Log marginal likelihood**: the value that a pass reports is not exact on a graph with loops
- **Parents at incompatible times**: when the messages from the parents of a node and its subtree leave no time with any probability, the node takes its time from its subtree alone, and the forward pass names these nodes in one warning. A reassortant whose segment trees disagree on its age gives such a node

## Depends on

- [M-optimize-segment-graphs-topology-fixed.md](M-optimize-segment-graphs-topology-fixed.md): the joint branch lengths that the time passes start from
