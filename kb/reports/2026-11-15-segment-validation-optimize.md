# Segment validation: branch-length optimization on simulated reassortment graphs

Validation of `optimize --segment-trees --segment-alignments` on segmented genomes simulated with reassortment, with [validation/segments.py](../../validation/segments.py) and [validation/reassortment.py](../../validation/reassortment.py) (README section in [validation/README.md](../../validation/README.md)).

## Setup

- **Graph**: coalescent with reassortment (CoalRe model; `Tc = 10` years, 80 tips at 40 sample times), reassortment rate 0.02 (`low`) or 0.2 (`high`) per lineage and year, 3 replicates each
- **Segments**: `s1`-`s4` with lengths 1500, 1200, 1000, 800 and rates 4, 3, 2, 1 x 10^-3 per site and year, each with a random GTR; sequences evolve along each segment's local tree
- **Input trees**: the true segment topologies with named internal nodes and branch lengths mu t x lognormal(0, 0.3); `full` keeps the nodes a segment passes through where another segment branches, `collapsed` drops them as tree inference would. Neither keeps the hidden nodes, which have one child in every segment (reassortment nodes, mergers of lineages with disjoint segments): no tree inference recovers them. Before commit `fcd860c10`, the `full` trees kept them
- **Variants**: `joint-full` and `joint-collapsed` (`optimize --segment-trees`), `separate` (`optimize --tree` per segment, topology kept)
- **Branches**: those of the collapsed segment trees; a fitted branch is the sum along its path in the output tree, times the fitted rate of the segment, times L. The two branches at each segment root are scored as a sum, because their split is not identifiable

## Results (commit `fcd860c10`)

Reassortments per replicate: high 27, high 16, high 27, low 0, low 3, low 5

## Relative segment rates (length-weighted mean 1)

|                           |   true |   estimated |    sd |
|:--------------------------|-------:|------------:|------:|
| ('joint-collapsed', 's1') |  1.452 |       1.461 | 0.04  |
| ('joint-collapsed', 's2') |  1.089 |       1.088 | 0.016 |
| ('joint-collapsed', 's3') |  0.726 |       0.703 | 0.043 |
| ('joint-collapsed', 's4') |  0.363 |       0.376 | 0.053 |
| ('joint-full', 's1')      |  1.452 |       1.462 | 0.041 |
| ('joint-full', 's2')      |  1.089 |       1.086 | 0.021 |
| ('joint-full', 's3')      |  0.726 |       0.705 | 0.045 |
| ('joint-full', 's4')      |  0.363 |       0.374 | 0.049 |

## Branch error against expected substitutions

RMSE in substitutions of the fitted substitutions (rate x length x L) against the expected substitutions mu t L, which a fit pooled over segments targets, over the branches of the collapsed segment trees, by the number of segments whose collapsed tree contains the same branch. `poisson` is the RMSE of Poisson counts around mu t L.

|             |   branches |   separate |   joint-collapsed |   joint-full |   poisson |
|:------------|-----------:|-----------:|------------------:|-------------:|----------:|
| ('high', 1) |        149 |       2.31 |              2.11 |         2.08 |      2.34 |
| ('high', 2) |        212 |       1.82 |              1.37 |         1.38 |      1.93 |
| ('high', 3) |        171 |       1.86 |              1.49 |         1.48 |      1.89 |
| ('high', 4) |       1340 |       1.23 |              0.7  |         0.7  |      1.17 |
| ('low', 1)  |          8 |       3.11 |              2.81 |         2.28 |      2.2  |
| ('low', 2)  |         28 |       2.35 |              1.93 |         1.93 |      1.97 |
| ('low', 3)  |         12 |       2.3  |              1.7  |         1.62 |      1.32 |
| ('low', 4)  |       1824 |       1.51 |              0.81 |         0.81 |      1.48 |

## Branch error against realized substitutions

RMSE in substitutions of the fitted substitutions (rate x length x L) against the segment's own realized substitutions, which a separate fit follows, over the branches of the collapsed segment trees, by the number of segments whose collapsed tree contains the same branch. `poisson` is the RMSE of Poisson counts around mu t L.

|             |   branches |   separate |   joint-collapsed |   joint-full |   poisson |
|:------------|-----------:|-----------:|------------------:|-------------:|----------:|
| ('high', 1) |        149 |       0.76 |              1.22 |         1.29 |      2.34 |
| ('high', 2) |        212 |       0.38 |              1.28 |         1.28 |      1.93 |
| ('high', 3) |        171 |       0.26 |              1.3  |         1.3  |      1.89 |
| ('high', 4) |       1340 |       0.1  |              0.97 |         0.97 |      1.17 |
| ('low', 1)  |          8 |       0.39 |              1.02 |         1.46 |      2.2  |
| ('low', 2)  |         28 |       0.66 |              1.33 |         1.33 |      1.97 |
| ('low', 3)  |         12 |       0.45 |              0.9  |         0.9  |      1.32 |
| ('low', 4)  |       1824 |       0.34 |              1.28 |         1.28 |      1.48 |

## Sum of the two branches at each segment root against expected substitutions (RMSE)

| case   |   separate |   joint-collapsed |   joint-full |
|:-------|-----------:|------------------:|-------------:|
| high   |      10.57 |             10.49 |        10.38 |
| low    |       8.24 |              4.8  |         5.1  |

## Findings

- The relative segment rates are recovered within a few percent
- Branches shared by several segments are fitted closer to the expected substitutions than by separate fits, the more so the more segments share them (4 segments: RMSE 0.7-0.8 against 1.2-1.5 substitutions, below the Poisson noise of one segment); branches of one segment are fitted about as well as by separate fits. Against a segment's own realized substitutions separate fits are closer by construction, since a shared length pools the segments
- The first runs on `full` trees exposed non-converging optimization (the log-likelihood fell by up to 219 units per iteration): the loop updated all edges at once, which overshoots on pairs of edges that a segment sees only through their sum (through a pass-through node, below a segment root) and along unidentifiable node times. The fix groups such edges ([kb/algo/optimization.md](../algo/optimization.md#branch-lengths-on-graphs-of-segment-trees)); after it (on the earlier `full` trees with hidden nodes), 11 of 12 runs increased the likelihood in every iteration, and one dipped by 4 units in the partition-rate step
- The hidden nodes caused the outlier at the segment roots of the earlier `full` trees (root-sum RMSE 22.2 substitutions in `high`, the root sum of `s1` in `high/rep02` 180 against 249 expected). Optimize can only hold the edge below such a node at length 0, which moved length between the root edges. Without hidden nodes `full` and `collapsed` agree (root-sum RMSE 10.4 and 10.5)
- On `full` trees, pass-through insertion adds 1 node in one replicate although the trees already hold every pass-through node: a node that one segment reaches directly and another through an intermediate node can be a reassortment, which topology alone does not distinguish ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name))
