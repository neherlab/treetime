# Hidden nodes of segment graphs are kept

A hidden node is a node with one child in every segment of a segment graph: a reassortment node (two parents and one child, each parent edge carrying part of the segments), a merger of two lineages that carry disjoint sets of segments, or one of the two nodes of a reassortment whose lineages merge again at one node. No segment branches there, so no tree inference recovers such a node, and the sequences carry no information about its date. TreeTime keeps hidden nodes as given: it neither contracts them nor rejects an input that has them.

**Type**: Scientific default of a new v1 input (documented limitation). Approved.

Validation: [kb/reports/2026-11-15-segment-validation-timetree.md](../reports/2026-11-15-segment-validation-timetree.md) (section Hidden nodes), with `validation/segments.py`.

## Effect on results

Simulated reassortment graphs (CoalRe, `Tc = 10` years, 80 tips over 4 years, four segments, reassortment rate 0.02 (`low`) or 0.2 (`high`) per lineage and year, 3 replicates each), dated by `timetree --segment-trees` at the true clock rate with `--confidence`. The same segment trees are given in three forms: `kept` with the hidden nodes, `contracted` without them (an edge across a hidden node gets the sum of the lengths of its segment), and `collapsed`, without hidden and pass-through nodes, as tree inference gives them. At the nodes that all three forms date, which leaves out the hidden nodes:

| case | kind            | form       | nodes | bias (y) | RMSE (y) | 90% coverage |
| ---- | --------------- | ---------- | ----: | -------: | -------: | -----------: |
| high | one parent      | kept       |   263 |   -0.184 |    0.894 |        0.829 |
| high | one parent      | contracted |   263 |   -0.104 |    0.573 |        0.863 |
| high | one parent      | collapsed  |   263 |   -0.093 |    0.565 |        0.878 |
| high | several parents | kept       |    32 |   -0.134 |    0.494 |        0.812 |
| high | several parents | contracted |    32 |   -0.135 |    0.575 |        0.812 |
| high | several parents | collapsed  |    32 |   -0.084 |    0.643 |        0.906 |
| low  | one parent      | kept       |   238 |   -0.125 |    0.839 |        0.861 |
| low  | one parent      | contracted |   238 |   -0.051 |    0.216 |        0.870 |
| low  | one parent      | collapsed  |   238 |   -0.053 |    0.203 |        0.870 |
| low  | several parents | kept       |     5 |   -0.178 |    0.250 |        0.800 |
| low  | several parents | contracted |     5 |   -0.190 |    0.255 |        0.800 |
| low  | several parents | collapsed  |     5 |   -0.189 |    0.258 |        0.800 |

- **Hidden nodes**: their 90% regions cover about half of the true dates (high: 0.54 of 69 reassortment nodes, 0.50 of 10 mergers; low: 0.50 of 6 reassortment nodes and of 4 loop nodes), and their RMSE is 1.5 to 5.9 years
- **Other nodes**: nodes without an edge to a hidden node move by 0.023 years RMS between `kept` and `contracted`. The harm is at the deep neighbours of hidden nodes, near the roots of the graph: errors up to 11.9 years in `kept` against 1.1 years in `contracted` (`low/rep01` `node_0085`, 36 years before the last sample), and shifts between the two forms up to 13.1 years
- **Several parents**: the `kept` row of `high` is lower than the `contracted` row, on 32 nodes; this is not a measured advantage of `kept`

## Mechanism

- **No data on the hidden node**: each segment passes the hidden node with one child, so the likelihood of a segment depends only on the sum of the edges along its path from a parent to the child
- **The edge below a hidden node is fitted at about zero**: the optimize step leaves 73 of 100 edges below hidden nodes shorter than 10^-8; the edges above them take the rest of the path. The input lengths are kept as the start (`InitialGuessMode::Auto`), so the zeros come from the optimizer. Why the optimizer ends at zero is unverified: no one traced it in the code
- **Paths through hidden nodes stretch**: the time passes give each edge the likelihood of its path as a function of its own length, with the other edges at their fitted lengths. An edge at length zero gets a narrow likelihood at zero, which places the hidden node at its child. On paths through one hidden node, the inferred time span exceeds the fitted path length by 0.27 years on average, against 0.10 years for the same parent and child after contraction; paths through two or three hidden nodes stretch by 1 to 17 years. That the cut-off of these likelihoods at zero causes the stretch is unverified: the numbers show the stretch, and no one traced its cause in the code

## Guidance to users

To keep the dates of the other nodes unaffected, contract the hidden nodes on each segment tree before the run: remove the node, and give its child the edge to its parent with the sum of the two lengths. Contraction on the input trees gives the same edges and scopes as contraction in the merged graph only before the merge inserts pass-through nodes ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)); after that step the two differ in 1 of the 6 simulated replicates, because the merge sees other paths when the hidden nodes are present. Inputs from tree inference and TreeKnit have no hidden nodes.

## Options not chosen

- **Contraction on each segment tree before the merge, with a warning that names the contracted nodes**: it gives the other nodes the dates of the input without hidden nodes and loses no information, but TreeTime would remove nodes that the user named, and the dates of these nodes would be missing from the outputs
- **Rejection of the input**: it changes no input without the user knowing, but every user with such an input would have to contract the nodes by hand to run TreeTime at all
