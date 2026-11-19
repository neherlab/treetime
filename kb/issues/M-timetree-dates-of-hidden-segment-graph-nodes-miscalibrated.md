# Dates of hidden nodes in segment graphs are miscalibrated, and paths through them stretch

A hidden node is a node with one child in every segment of a segment graph. TreeTime keeps such nodes as given ([kb/decisions/segment-graph-hidden-nodes-kept.md](../decisions/segment-graph-hidden-nodes-kept.md)). The sequences carry no information about their dates: every segment passes the node with one child, so the likelihood of a segment depends only on the sum of the edges along its path from a parent to the child. The time passes nevertheless give these nodes narrow regions, and they lengthen the paths through them.

Measured on simulated reassortment graphs ([kb/reports/2026-11-15-segment-validation-timetree.md](../reports/2026-11-15-segment-validation-timetree.md), section "Hidden nodes", `kept` form):

- **Hidden nodes**: the 90% regions cover about half of the true dates (high: 0.54 of 69 reassortment nodes, 0.50 of 10 mergers; low: 0.50 of 6 reassortment nodes and of 4 loop nodes), with RMSE 1.5 to 5.9 years and a bias of +0.43 years towards the child at reassortment nodes
- **Paths**: through one hidden node, the inferred time span exceeds the fitted path length by 0.27 years on average, against 0.10 years for the same parent and child after contraction. Through two or three hidden nodes, paths stretch by 1 to 17 years
- **Deep neighbours**: errors up to 11.9 years, against 1.1 years with the hidden nodes contracted (`low/rep01` `node_0085`); RMSE at nodes with one parent 0.894 and 0.839 years, against 0.573 and 0.216 years contracted

The optimize step leaves 73 of 100 edges below hidden nodes shorter than 10^-8, and the edges above them take the rest of the path. The time passes then give each edge the likelihood of its path as a function of its own length, with the other edges at their fitted lengths, so an edge at length zero gets a narrow likelihood at zero, which places the hidden node at its child.

> [!IMPORTANT]
> **Investigation required.** Why the optimizer ends these edges at zero: the input lengths are the start (`InitialGuessMode::Auto`), so the zeros come from the optimizer, and no one traced the path in the code. Evidence to collect: trace the fit of the two edges of one path through a hidden node of `low/rep01`.

> [!IMPORTANT]
> **Investigation required.** Whether the cut-off of the edge likelihoods at zero causes the stretch: the numbers show the stretch, and no one traced its cause in the time passes. Evidence to collect: the branch-length likelihoods that the time passes build for the edges of the chain `node_0085 -> node_0084 -> node_0083 -> node_0082` of `low/rep01` (14.3 years of time span for 5.4 years of fitted length), against the likelihood of the summed path.

The same unidentified split occurs below a segment root, where the refinement rounds reduce its effect, and in the series chains of [M-optimize-segment-graphs-topology-fixed.md](M-optimize-segment-graphs-topology-fixed.md). A fix that dates the edges through a node without data by the likelihood of their summed length changes scientific output and needs an entry in `kb/decisions/`.

## Validation

`validation/segments.py` on the `kept` form: the 90% regions of hidden nodes cover about 90% of the true dates, and the RMSE at the other nodes approaches that of the `contracted` form.
