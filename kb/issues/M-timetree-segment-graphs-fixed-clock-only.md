# Timetree on graphs of segment trees runs only at a fixed clock rate

`timetree --segment-trees --segment-alignments --clock-rate` dates a graph of segment trees merged by node name ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)). It optimizes the branch lengths jointly over the segments, then runs the time passes with several parents per node as loopy message passing, with refinement rounds that feed the clock branch lengths back into the marginal reconstruction ([kb/algo/timetree.md](../algo/timetree.md#graphs-with-several-parents)). Validation: [kb/reports/2026-11-15-segment-validation-timetree.md](../reports/2026-11-15-segment-validation-timetree.md).

`fn run_segment_timetree` [packages/treetime/src/timetree/segments.rs](../../packages/treetime/src/timetree/segments.rs) is a reduced pipeline. Compared to the tree pipeline it lacks:

- **Clock rate estimation**: `--clock-rate` is required, and `--clock-std-dev` is rejected. The root-to-tip regression needs one divergence per tip, which a tip with several root paths does not have; the marginal likelihood probe of the rate should carry over
- **Coalescent prior**: rejected. The expected lineage count of a segment graph is per segment, and the prior's merger events are the branching events of each segment
- **Rooting**: the roots of the graph are kept; `--reroot` and `--outgroup` are rejected. Rerooting inverts edges, which has no meaning at a node with several parents
- **Polytomy resolution, short-branch collapse, clock filter, tip pruning, relaxed clock**: not run, or rejected
- **Outputs**: reconstructed sequences, augur node data, the clock model and the root-to-tip table are not written; the node dates table covers every node of the graph

Known limitations of the passes:

- **Hidden nodes**: a node with one child in every segment (a reassortment node given as a pass-through node of every segment, or a merger of lineages with disjoint segments) is placed at its child with a narrow region, because optimize holds the edge below it at length 0. Contracting such nodes for the time inference, or reporting them without a date, would be honest
- **Zero-length edges**: edges of length 0 from different polytomy resolutions in the segment trees close spurious loops (the flu graph has 391 nodes with a zero-length parent edge among 801 with several parents). Contracting them needs the scopes to follow the contraction
- **Sweeps**: the sweep limit of 20 is reached on the flu graph at a relative change of 2 x 10^-4; the log marginal likelihood that a sweep reports is not exact on a graph with loops
- **Default of the cross-parent messages**: exact on small loops, no measurable effect on the simulated graphs and the flu graph, twice the run time; `--no-cross-parent-messages` switches them off

## Depends on

- [M-optimize-segment-graphs-topology-fixed.md](M-optimize-segment-graphs-topology-fixed.md): the joint branch lengths that the time passes start from
