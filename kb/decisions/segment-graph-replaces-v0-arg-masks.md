# Segment trees run on one graph instead of v0's focal tree with masks

v0 `treetime arg` dates one focal segment tree at a time. The other segments count only on the branches inside a maximum compatible clade (MCC) of TreeKnit and are masked elsewhere. v1 merges all segment trees into one graph by node name and records on each edge which segments it carries (edge scopes). Every segment then runs on its own complete tree, every shared edge has one branch length, and every node has one date. Results of `ancestral`, `optimize` and `timetree` with `--segment-trees` therefore differ from v0 `arg` by design, and v0 is not a test oracle for segment graphs.

**Type**: Intentional scientific divergence from v0 (new v1 model of reassortment). Approved.

## v0 masks

v0 gives every branch of the focal tree a mask: the whole alignment when parent and child lie in the same MCC, the focal segment only otherwise (`fn setup_arg()` [packages/legacy/treetime/treetime/arg.py#L128-L133](../../packages/legacy/treetime/treetime/arg.py#L128-L133)). The marginal backward pass sends a uniform message across a branch at masked sites ([packages/legacy/treetime/treetime/treeanc.py#L862-L871](../../packages/legacy/treetime/treetime/treeanc.py#L862-L871)), the forward pass replaces the parent message at masked sites by one ([packages/legacy/treetime/treetime/treeanc.py#L909-L916](../../packages/legacy/treetime/treetime/treeanc.py#L909-L916)), and the masks enter branch-length optimization and the timetree through the pattern multiplicity ([packages/legacy/treetime/treetime/clock_tree.py#L363](../../packages/legacy/treetime/treetime/clock_tree.py#L363)).

Consequences in v0:

- **Cut local trees**: the tree of every non-focal segment is cut at each MCC root, and the cut subtrees have no root prior
- **One answer per focal tree**: a node that several segments share gets a date from each focal-tree run, and the runs need not agree
- **Two trees**: the command takes exactly two trees and one MCC file

## v1 graph

- **Merge by node name**: a name that occurs in several trees is one graph node, a parent-child pair that occurs in several trees is one edge, and the merge inserts the nodes that a segment only passes through ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name))
- **Edge scopes**: each edge carries a set of segments, and each site has at most one parent at every node ([packages/treetime-graph/src/scope.rs](../../packages/treetime-graph/src/scope.rs))
- **Exact reconstruction per segment**: the Fitch and marginal passes send each site only along the edges that carry it, so every site sees its own tree, and each segment gives the same sequences, mutations and log-likelihood as a run on its own tree, also with gamma rate categories and with several partitions per segment ([packages/treetime/src/ancestral/__tests__/test_pipeline_segment_graph.rs](../../packages/treetime/src/ancestral/__tests__/test_pipeline_segment_graph.rs))
- **Shared parameters**: an edge that several segments carry has one length, fitted to all of them, and a node has one date for all segments ([timetree-segment-graph-loopy-message-passing.md](timetree-segment-graph-loopy-message-passing.md))

## Consequences

- **Oracles**: tests compare each segment with a run on its own tree, and dates with the truth of simulated graphs ([kb/reports/2026-11-15-segment-validation-timetree.md](../reports/2026-11-15-segment-validation-timetree.md)), not with v0
- **Input**: named segment trees, one per segment, in place of two trees and an MCC file. Turning TreeKnit output into named trees is open ([kb/reports/2026-11-15-treeknit-output-to-named-segment-trees.md](../reports/2026-11-15-treeknit-output-to-named-segment-trees.md))
- **Command**: the v0 command line `treetime arg` stops with a message that names `treetime timetree --segment-trees --segment-alignments --clock-rate` ([packages/app-cli/src/cli/v0_args.rs](../../packages/app-cli/src/cli/v0_args.rs))

Design and alternatives: [kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md).
