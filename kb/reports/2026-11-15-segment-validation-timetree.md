# Segment validation: dates on simulated reassortment graphs

Validation of `timetree --segment-trees --segment-alignments --clock-rate` on the simulated segment graphs of [kb/reports/2026-11-15-segment-validation-optimize.md](2026-11-15-segment-validation-optimize.md) (same data, same seeds), with [validation/segments.py](../../validation/segments.py).

## Setup

- **Data**: as in the optimize report: CoalRe graphs with `Tc = 10` years, 80 tips sampled over 4 years, reassortment rate 0.02 (`low`) or 0.2 (`high`) per lineage and year, 3 replicates each; four segments with rates 4, 3, 2, 1 x 10^-3 per site and year
- **Dates**: exact tip dates from the sampling schedule; fixed clock rate 2.756 x 10^-3, the true length-weighted mean rate; 90% maximum posterior regions (`--confidence`); the default 2 refinement rounds
- **Variants**: the `full` trees (with the nodes a segment passes through where another segment branches) and the `collapsed` trees (without them, as tree inference gives; the merge inserts the pass-through nodes it can determine), each without cross-parent messages (the default) and with `--cross-parent-messages` (`cross`)
- **Hidden nodes**: neither form contains nodes with one child in every segment (see Findings); before commit `f16b1235a` the `full` trees kept them
- **Node kinds**: `several parents` counts nodes with more than one parent in the graph of the variant; the collapsed trees, which also drop the pass-through nodes, give more of them

## Results (commit `f16b1235a`)

| case | variant          | kind            | nodes | bias (y) | RMSE (y) | coverage | median width (y) | sweeps |
| ---- | ---------------- | --------------- | ----: | -------: | -------: | -------: | ---------------: | -----: |
| high | collapsed        | one parent      |   229 |   -0.083 |    0.638 |    0.891 |            0.395 |      1 |
| high | collapsed        | several parents |    66 |   -0.120 |    0.246 |    0.848 |            0.493 |      1 |
| high | collapsed, cross | one parent      |   229 |   -0.078 |    0.637 |    0.891 |            0.395 |    4.0 |
| high | collapsed, cross | several parents |    66 |   -0.116 |    0.241 |    0.848 |            0.493 |    4.0 |
| high | full             | one parent      |   263 |   -0.104 |    0.573 |    0.863 |            0.398 |      1 |
| high | full             | several parents |    32 |   -0.135 |    0.575 |    0.812 |            0.583 |      1 |
| high | full, cross      | one parent      |   263 |   -0.102 |    0.575 |    0.867 |            0.398 |    3.7 |
| high | full, cross      | several parents |    32 |   -0.117 |    0.581 |    0.812 |            0.582 |    3.8 |
| low  | collapsed        | one parent      |   232 |   -0.052 |    0.203 |    0.866 |            0.398 |      1 |
| low  | collapsed        | several parents |    11 |   -0.138 |    0.221 |    0.909 |            0.760 |      1 |
| low  | collapsed, cross | one parent      |   232 |   -0.052 |    0.203 |    0.866 |            0.398 |    2.6 |
| low  | collapsed, cross | several parents |    11 |   -0.135 |    0.224 |    0.909 |            0.759 |    3.8 |
| low  | full             | one parent      |   238 |   -0.051 |    0.216 |    0.870 |            0.401 |      1 |
| low  | full             | several parents |     5 |   -0.190 |    0.255 |    0.800 |            0.756 |      1 |
| low  | full, cross      | one parent      |   238 |   -0.051 |    0.216 |    0.870 |            0.401 |    2.3 |
| low  | full, cross      | several parents |     5 |   -0.189 |    0.255 |    0.800 |            0.752 |    3.0 |

Sweeps are those of the first round; without cross-parent messages one sweep is exact for the passes, because the backward pass does not read the previous sweep.

## Findings

- **Collapsed trees, the realistic input, are dated well**: bias below 0.14 years, RMSE 0.2 to 0.64 years, and the 90% regions cover 85 to 91% of the true dates, at nodes with one parent and at nodes with several parents alike
- **Cross-parent messages make no measurable difference here**: on both forms the two variants agree to within 0.01 years in RMSE. On the exact small-loop tests ([`packages/treetime/src/timetree/inference/__tests__/test_reassortment_passes.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_reassortment_passes.rs)) they differ by 10x in CDF error, so the loops of these graphs must be long or weakly informative. On the flu graph (below) the dates of the two variants differ by less than 0.09 years at 99% of nodes
- **Decision**: cross-parent messages are off by default. They carry little information in practice and slow convergence (up to 20 sweeps per round on flu); `--cross-parent-messages` switches them on
- **Full and collapsed trees agree**: bias below 0.2 years, RMSE 0.2 to 0.64 years, coverage 0.80 to 0.91, at nodes with one parent and with several parents
- **Hidden nodes**: nodes with one child in every segment leave no branching in any segment, so no tree inference recovers them and sequences do not date them. In the simulated graphs (high, 3 replicates) they are 69 reassortment nodes (two parents, one child, each parent edge carrying part of the segments), 10 mergers of lineages that carry disjoint sets of segments, and 2 nodes of one reassortment whose two lineages merge again at the same node; low: 6, 0 and 4. The CoalRe model allows all three. Dated as nodes of the graph, they sit near their child with a narrow region and harm the dates of their deep neighbours; contracting them before the merge removes the harm (section Hidden nodes below)
- **Refinement rounds are needed**: optimize fixes only the sum of the two edges below a segment root. With a single time pass (before commit `89c301db7`, on the trees with hidden nodes), the arbitrary split moved deep nodes by up to 21 years; the refinement rounds feed the clock branch lengths back into the marginal reconstruction. The largest remaining errors, 6.3 years (collapsed) and 7.1 years (full) in `high`, sit at the oldest nodes, near the roots of the graph

## Hidden nodes

Comparison of three inputs with the same simulated data (commit `f16b1235a` changed no random draw): `kept`, the `full` trees with the hidden nodes as nodes of one child in every segment; `contracted`, the `full` trees without them, where an edge across a hidden node gets the sum of the perturbed lengths of its segment (the `timetree-full` variant); and `collapsed`. Fixed clock rate, `--confidence`, default message passing, 3 replicates per case. `validation/segments.py` writes the `kept` trees and runs `timetree-kept` and `joint-kept`.

TreeTime keeps hidden nodes as given, based on these measurements: [kb/decisions/segment-graph-hidden-nodes-kept.md](../decisions/segment-graph-hidden-nodes-kept.md).

The runs use the binary of commit `c06d8a04f`, which reproduces the `full` rows of the results above to the last digit. Commit `c4c2b5e9a` changed the merge after it, and a later fix restored its rule (section Merge of trees that hold all their pass-through nodes). With the fix, the code gives the tables below to the last digit for `kept` and `contracted`, and the same `collapsed` dates as commit `873490a98`; `collapsed` differs from `c06d8a04f` by at most 0.002 years in bias and RMSE, through other commits.

### Results (commit `c06d8a04f`)

Shared nodes are the nodes that all three graphs date, which leaves out the hidden nodes. The kind of a node is its kind in the contracted graph, so each case and kind has the same nodes in every form.

| case | kind            | form       | nodes | bias (y) | RMSE (y) | coverage | median width (y) |
| ---- | --------------- | ---------- | ----: | -------: | -------: | -------: | ---------------: |
| high | one parent      | kept       |   263 |   -0.184 |    0.894 |    0.829 |            0.412 |
| high | one parent      | contracted |   263 |   -0.104 |    0.573 |    0.863 |            0.398 |
| high | one parent      | collapsed  |   263 |   -0.093 |    0.565 |    0.878 |            0.399 |
| high | several parents | kept       |    32 |   -0.134 |    0.494 |    0.812 |            0.593 |
| high | several parents | contracted |    32 |   -0.135 |    0.575 |    0.812 |            0.583 |
| high | several parents | collapsed  |    32 |   -0.084 |    0.643 |    0.906 |            0.608 |
| low  | one parent      | kept       |   238 |   -0.125 |    0.839 |    0.861 |            0.403 |
| low  | one parent      | contracted |   238 |   -0.051 |    0.216 |    0.870 |            0.401 |
| low  | one parent      | collapsed  |   238 |   -0.053 |    0.203 |    0.870 |            0.401 |
| low  | several parents | kept       |     5 |   -0.178 |    0.250 |    0.800 |            0.755 |
| low  | several parents | contracted |     5 |   -0.190 |    0.255 |    0.800 |            0.756 |
| low  | several parents | collapsed  |     5 |   -0.189 |    0.258 |    0.800 |            0.759 |

Hidden nodes in `kept`; a loop node is one of the two nodes of a reassortment whose lineages merge again at one node:

| case | kind            | nodes | bias (y) | RMSE (y) | coverage | median width (y) |
| ---- | --------------- | ----: | -------: | -------: | -------: | ---------------: |
| high | reassortment    |    69 |    0.425 |    1.517 |    0.536 |            0.759 |
| high | disjoint merger |    10 |   -1.478 |    2.824 |    0.500 |            2.714 |
| high | loop            |     2 |    0.042 |    0.110 |    1.000 |            0.429 |
| low  | reassortment    |     6 |    1.860 |    3.212 |    0.500 |            0.954 |
| low  | loop            |     4 |   -4.723 |    5.904 |    0.500 |            8.346 |

Date shift `kept` minus `contracted` at shared nodes. A neighbour is linked to a hidden node by an edge or through other hidden nodes; age in years before 2025.0:

| case | neighbour | age (y) | nodes | RMS shift (y) | largest shift (y) | RMSE kept (y) | RMSE contracted (y) |
| ---- | --------- | ------- | ----: | ------------: | ----------------: | ------------: | ------------------: |
| high | no        | 0-10    |   170 |         0.023 |              0.11 |         0.151 |               0.142 |
| high | yes       | 0-5     |    80 |         0.099 |              0.39 |         0.224 |               0.200 |
| high | yes       | 5-10    |    29 |         0.281 |              0.68 |         0.677 |               0.672 |
| high | yes       | 10-20   |     7 |         0.956 |              1.47 |         0.942 |               1.243 |
| high | yes       | 20+     |     9 |         3.685 |              6.43 |         4.605 |               2.717 |
| low  | no        | 0-20    |   228 |         0.023 |              0.33 |         0.159 |               0.154 |
| low  | yes       | 0-20    |    12 |         0.060 |              0.14 |         0.603 |               0.591 |
| low  | yes       | 20+     |     3 |         7.888 |             13.07 |         7.257 |               0.785 |

The largest shifts: `low/rep01` `node_0085` (age 36 years) -13.1 years, error -11.9 kept and +1.1 contracted; `high/rep00` `node_0133` and `node_0132` (age 25 to 26 years) -6.4 and -5.7 years; `high/rep02` `node_0131` (age 31 years) +4.3 years. All lie above chains of hidden nodes near the roots of the graph.

### Mechanism

- **The sequences carry no information about a hidden node**: every segment passes it with one child, so the likelihood of a segment depends only on the sum of the edges along its path from a parent to the child. Moving the hidden node changes no likelihood; its date is bounded only by its child and its parents
- **Optimize ends most edges below a hidden node at length zero**: the joint fit (`joint-kept`) gives 84% (high) and 67% (low) of the per-segment edges below a hidden node less than 0.01 substitutions, median 0, where the expected median is 0.93 and 6.1 substitutions; the edges above them take the rest of the path (median 2.6 and 8.0 substitutions against 1.7 and 1.9 expected). The optimize step of `timetree-kept` leaves 73 of 100 such edges below 10^-8. The split is not fixed by the data, so other edges of the same kind end elsewhere: in `low/rep01` the loop chain `node_0085 -> node_0084 -> node_0083 -> node_0082` gets 0.0077, 0.0061 and 0.0011
- **The time passes treat the edges of a path as independent branch-length likelihoods**: each edge gets the likelihood of its path as a function of its own length, with the other edges at their fitted lengths. An edge at length zero gets a narrow likelihood at zero, which places the hidden node at its child (reassortment nodes: bias +0.43 years towards the child, coverage 0.54, median width 0.76 years). Each such edge adds positive length to the path, because its likelihood is cut off at zero: on paths with one hidden node the inferred time span exceeds the fitted path length by 0.27 years (median 0.16) in `kept` against 0.10 years (median 0.05) for the same parent and child in `contracted`, and the error against the true span is +0.22 against -0.03 years. Paths with two or three hidden nodes stretch by 1 to 17 years: the chain of `low/rep01` above spans 14.3 years against 5.4 years of fitted length and 6.1 true years
- **Contraction gives the constraint that the data give**: after contraction, the child has the parent edges of the hidden node, each with the sum of its path. The data fix the same sums; nothing fixes the hidden node, and contraction removes it with no information lost. The deep neighbours then have the errors of the collapsed input

### Contraction of the input and contraction in the graph

- **Same topology and scopes before the pass-through step**: contracting each hidden node in the graph of the `kept` trees, where every child takes over the parent edges whose segments it shares, gives the edges and scopes of the `contracted` trees in all 6 replicates. The lengths differ: the input sums the lengths of each segment and the merge then takes the mean over segments, while contraction in the merged graph would sum the means. Optimize refits the lengths, so this changes only its start
- **Different after the pass-through step**: the merge (`insert_pass_through_nodes()`) sees other paths when the hidden nodes are present. After the merge, the contracted `kept` graph and the `contracted` graph differ by one edge in 1 of 6 replicates (`high/rep00`). So a contraction inside TreeTime belongs before the merge, on each segment tree

### Merge of trees that hold all their pass-through nodes

The `full` and `kept` trees are the true local trees, with every node that a segment passes through where another segment branches, so the merge should change no edge. It reroutes one edge in each form over the 6 replicates (`node_0132 -> node_0126` in `full` of `high/rep00`, `node_0017 -> node_0003` in `kept` of `high/rep02`): a segment edge `X -> Y` moves below a node `W` that the other segments put on the path and that the segment does not hold, which deletes a true parent edge of `Y`. Topology alone does not show that a reassortment separates the segment from `W` there ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)). `validation/segments.py evaluate` lists these edges in the table Merged graphs of `summary.md`.

Commit `c4c2b5e9a` also let the edge continue through nodes with one child in the segment's own input, until a fix limited the chain to nodes that the merge inserted or attached a child to. A node with one child in the input is a node that the segment passes through towards another node, and the direct edge `X -> Y` states that the segment bypasses it. With the chain through such nodes, the merge rerouted 5 edges in each form (`high/rep01`, `high/rep02`, `low/rep02`), and:

- `high/rep00` `contracted` failed: "Time inference needs a grid of 1804164 points, more than the limit of 1000000 ... when sending the time message forward to node 'node_0132'", with parent time distributions over [177.6, 2021.1] years
- `contracted` lost accuracy: RMSE at nodes with one parent 0.573 to 1.140 years (high, `rep01` and `rep02` only) and 0.216 to 0.317 years (low); `high/rep02` `node_0131` moved from -0.31 to -11.6 years, and in `low/rep02`, where the merge rerouted the s3 edge `node_0086 -> node_0080`, the parent `node_0088` of `node_0086` moved from +1.8 to +3.6 years
- `collapsed` was unchanged to the last digit

## Flu (four segments, 2839 nodes)

`timetree --segment-trees` on the named TreeKnit trees of [kb/reports/2026-11-15-treeknit-output-to-named-segment-trees.md](2026-11-15-treeknit-output-to-named-segment-trees.md) at `--clock-rate 0.003` with `--confidence`: 43 s and 350 MB with `--cross-parent-messages` (optimize about 15 s; 3 time rounds with 9, 20 and 16 sweeps; the second reached the sweep limit at a relative change of 2 x 10^-4), 19 s and 310 MB without (the default).

## Open

- **Merge of trees with pass-through nodes**: the merge moves an edge of a segment below a node that the other segments put on its path and that the segment does not hold, also where a reassortment separates them; in the simulation this changes one edge per form over 6 replicates (section Hidden nodes)
- **Sweep tolerance** (with `--cross-parent-messages`): the flu graph needs more than 20 sweeps in one round to reach 10^-4 of the date span
- Clock rate estimation, coalescent prior, rooting and polytomy resolution on segment graphs: [kb/issues/M-timetree-segment-graphs-fixed-clock-only.md](../issues/M-timetree-segment-graphs-fixed-clock-only.md)
