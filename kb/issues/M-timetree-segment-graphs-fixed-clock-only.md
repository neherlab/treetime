# Timetree on graphs of segment trees runs only at a fixed clock rate

`timetree --segment-trees --segment-alignments --clock-rate` dates a graph of segment trees merged by node name ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)). It optimizes the branch lengths jointly over the segments, then runs the time passes with several parents per node as loopy message passing, with refinement rounds that feed the clock branch lengths back into the marginal reconstruction ([kb/algo/timetree.md](../algo/timetree.md#graphs-with-several-parents)). Validation: [kb/reports/2026-11-15-segment-validation-timetree.md](../reports/2026-11-15-segment-validation-timetree.md).

`fn run_segment_timetree` [packages/treetime/src/timetree/segments.rs](../../packages/treetime/src/timetree/segments.rs) is a reduced pipeline. Compared to the tree pipeline it lacks:

- **Clock rate estimation**: `--clock-rate` is required, and `--clock-std-dev` is rejected. The root-to-tip regression needs one divergence per tip, which a tip with several root paths does not have; the marginal likelihood probe of the rate should carry over
- **Coalescent prior**: rejected. The expected lineage count of a segment graph is per segment, and the prior's merger events are the branching events of each segment
- **Rooting**: the roots of the graph are kept; `--reroot` and `--outgroup` are rejected. Rerooting inverts edges, which has no meaning at a node with several parents
- **Polytomy resolution, short-branch collapse, clock filter, tip pruning, relaxed clock**: not run, or rejected
- **Outputs**: reconstructed sequences, augur node data, the clock model and the root-to-tip table are not written; the node dates table covers every node of the graph
- **Branch-length settings**: the optimize step runs Brent with the fixed values `OPTIMIZE_MAX_ITER = 10` and `OPTIMIZE_TOLERANCE = 0.1` and no damping (`fn segment_optimize_params()` [packages/treetime/src/timetree/segments.rs#L201](../../packages/treetime/src/timetree/segments.rs#L201))

Each feature of the tree pipeline (`packages/treetime/src/timetree/pipeline.rs`) needs a second implementation in `run_segment_timetree`, or the two pipelines drift apart.

> [!IMPORTANT]
> **Decision required.** Extend `run_segment_timetree` feature by feature, or run segment graphs through the tree pipeline, with the steps that have no meaning on a graph with several parents (rerooting) rejected there.

Known limitations of the passes:

- **Hidden nodes**: kept as given, neither contracted nor rejected ([kb/decisions/segment-graph-hidden-nodes-kept.md](../decisions/segment-graph-hidden-nodes-kept.md)). A node with one child in every segment (a reassortment node, a merger of lineages with disjoint segments, or a node of a reassortment whose lineages merge again) is placed near its child with a narrow region (coverage about 0.5), and paths through it stretch, which moves its deep neighbours by up to 13 years. On simulated graphs the RMSE at nodes with one parent is 0.89 and 0.84 years with the hidden nodes, against 0.57 and 0.22 years with them contracted on each segment tree before the run. Tree inference and TreeKnit never produce such nodes. The miscalibration and its unverified causes: [M-timetree-dates-of-hidden-segment-graph-nodes-miscalibrated.md](M-timetree-dates-of-hidden-segment-graph-nodes-miscalibrated.md)
- **Zero-length edges**: kept, not contracted ([kb/decisions/segment-graph-zero-length-edges-kept.md](../decisions/segment-graph-zero-length-edges-kept.md)); inputs to TreeKnit should be collapsed trees
- **Cross-parent messages**: off by default ([kb/decisions/timetree-segment-graph-loopy-message-passing.md](../decisions/timetree-segment-graph-loopy-message-passing.md)); `--cross-parent-messages` switches them on. With them, the sweep limit of 20 is reached on the flu graph at a relative change of 2 x 10^-4
- **Log marginal likelihood**: the value that a pass reports is not exact on a graph with loops
- **Parents at incompatible times**: when the messages from the parents of a node and its subtree leave no time with any probability, the node takes its time from its subtree alone, and the forward pass names these nodes in one warning. A reassortant whose segment trees disagree on its age gives such a node

## Depends on

- [M-optimize-segment-graphs-topology-fixed.md](M-optimize-segment-graphs-topology-fixed.md): the joint branch lengths that the time passes start from
