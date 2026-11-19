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
- **Cross-parent messages make no measurable difference here**: on both forms the two variants agree to within 0.01 years in RMSE. On the exact small-loop tests ([packages/treetime/src/timetree/inference/__tests__/test_reassortment_passes.rs](../../packages/treetime/src/timetree/inference/__tests__/test_reassortment_passes.rs)) they differ by 10x in CDF error, so the loops of these graphs must be long or weakly informative. On the flu graph (below) the dates of the two variants differ by less than 0.09 years at 99% of nodes
- **Decision**: cross-parent messages are off by default. They carry little information in practice and slow convergence (up to 20 sweeps per round on flu); `--cross-parent-messages` switches them on
- **Full and collapsed trees agree**: bias below 0.2 years, RMSE 0.2 to 0.64 years, coverage 0.80 to 0.91, at nodes with one parent and with several parents
- **Hidden nodes**: nodes with one child in every segment leave no branching in any segment, so no tree inference recovers them and sequences do not date them. In the simulated graphs (high, 3 replicates) they are 69 reassortment nodes (two parents, one child, each parent edge carrying part of the segments), 10 mergers of lineages that carry disjoint sets of segments, and 2 reassortments whose two lineages merge again at the same node; low: 6, 0 and 4. The CoalRe model allows all three. With them in the `full` trees (before commit `f16b1235a`), they were placed at their child with a narrow region (coverage 0.5 to 0.6), because optimize holds the edge below them at length 0, and they pulled the neighbouring deep nodes (RMSE 0.6 to 1.3 years at the other nodes)
- **Refinement rounds are needed**: optimize fixes only the sum of the two edges below a segment root. With a single time pass (before commit `89c301db7`, on the trees with hidden nodes), the arbitrary split moved deep nodes by up to 21 years; the refinement rounds feed the clock branch lengths back into the marginal reconstruction. The largest remaining errors, 6.3 years (collapsed) and 7.1 years (full) in `high`, sit at the oldest nodes, near the roots of the graph

## Flu (four segments, 2839 nodes)

`timetree --segment-trees` on the named TreeKnit trees of [kb/reports/2026-11-15-treeknit-output-to-named-segment-trees.md](2026-11-15-treeknit-output-to-named-segment-trees.md) at `--clock-rate 0.003` with `--confidence`: 43 s and 350 MB with `--cross-parent-messages` (optimize about 15 s; 3 time rounds with 9, 20 and 16 sweeps; the second reached the sweep limit at a relative change of 2 x 10^-4), 19 s and 310 MB without (the default).

## Open

- **Sweep tolerance** (with `--cross-parent-messages`): the flu graph needs more than 20 sweeps in one round to reach 10^-4 of the date span
- Clock rate estimation, coalescent prior, rooting and polytomy resolution on segment graphs: [kb/issues/M-timetree-segment-graphs-fixed-clock-only.md](../issues/M-timetree-segment-graphs-fixed-clock-only.md)
