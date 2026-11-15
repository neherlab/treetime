# Segment validation: dates on simulated reassortment graphs

Validation of `timetree --segment-trees --segment-alignments --clock-rate` on the simulated segment graphs of [kb/reports/2026-11-15-segment-validation-optimize.md](2026-11-15-segment-validation-optimize.md) (same data, same seeds), with [validation/segments.py](../../validation/segments.py).

## Setup

- **Data**: as in the optimize report: CoalRe graphs with `Tc = 10` years, 80 tips sampled over 4 years, reassortment rate 0.02 (`low`) or 0.2 (`high`) per lineage and year, 3 replicates each; four segments with rates 4, 3, 2, 1 x 10^-3 per site and year
- **Dates**: exact tip dates from the sampling schedule; fixed clock rate 2.756 x 10^-3, the true length-weighted mean rate; 90% maximum posterior regions (`--confidence`); the default 2 refinement rounds
- **Variants**: the `full` trees (with the nodes a segment passes through) and the `collapsed` trees (without them, as tree inference gives; the merge inserts the pass-through nodes it can determine), each with cross-parent messages (default) and with `--no-cross-parent-messages`
- **Node kinds**: `hidden` nodes have one child in every segment: reassortment nodes, and mergers of lineages that carry disjoint sets of segments. No segment branches there, so sequences do not date them; the collapsed trees leave them out. `several parents` counts nodes with more than one parent in the graph of the variant

## Results (commit `09294169f`)

| case | variant                     | kind            | nodes | bias (y) | RMSE (y) | coverage | median width (y) | sweeps |
| ---- | --------------------------- | --------------- | ----: | -------: | -------: | -------: | ---------------: | -----: |
| high | collapsed                   | one parent      |   229 |   -0.078 |    0.637 |    0.891 |            0.395 |    4.0 |
| high | collapsed                   | several parents |    66 |   -0.116 |    0.241 |    0.848 |            0.493 |    4.0 |
| high | collapsed, children only    | one parent      |   229 |   -0.083 |    0.638 |    0.891 |            0.395 |      1 |
| high | collapsed, children only    | several parents |    66 |   -0.120 |    0.246 |    0.848 |            0.493 |      1 |
| high | full                        | hidden          |    81 |   -0.141 |    2.085 |    0.494 |            0.817 |    8.3 |
| high | full                        | one parent      |   295 |   -0.266 |    1.271 |    0.837 |            0.417 |    7.7 |
| high | full, children only         | hidden          |    81 |    0.180 |    1.716 |    0.543 |            0.820 |      1 |
| high | full, children only         | one parent      |   295 |   -0.179 |    0.860 |    0.827 |            0.437 |      1 |
| low  | collapsed                   | one parent      |   232 |   -0.052 |    0.203 |    0.866 |            0.398 |    2.6 |
| low  | collapsed                   | several parents |    11 |   -0.135 |    0.224 |    0.909 |            0.759 |    3.8 |
| low  | collapsed, children only    | one parent      |   232 |   -0.052 |    0.203 |    0.866 |            0.398 |      1 |
| low  | collapsed, children only    | several parents |    11 |   -0.138 |    0.221 |    0.909 |            0.760 |      1 |
| low  | full                        | hidden          |    10 |   -0.187 |    3.754 |    0.600 |            4.001 |    3.5 |
| low  | full                        | one parent      |   243 |   -0.107 |    0.597 |    0.860 |            0.406 |    2.7 |
| low  | full, children only         | hidden          |    10 |   -0.773 |    4.487 |    0.500 |            4.002 |      1 |
| low  | full, children only         | one parent      |   243 |   -0.126 |    0.831 |    0.860 |            0.406 |      1 |

Sweeps are those of the first round; with `--no-cross-parent-messages` one sweep is exact for the passes, because the backward pass does not read the previous sweep.

## Findings

- **Collapsed trees, the realistic input, are dated well**: bias below 0.14 years, RMSE 0.2 to 0.64 years, and the 90% regions cover 85 to 91% of the true dates, at nodes with one parent and at nodes with several parents alike
- **Cross-parent messages make no measurable difference here**: on the collapsed graphs the two variants agree to within 0.01 years in RMSE. On the exact small-loop tests ([packages/treetime/src/timetree/inference/__tests__/test_reassortment_passes.rs](../../packages/treetime/src/timetree/inference/__tests__/test_reassortment_passes.rs)) they differ by 10x in CDF error, so the loops of these graphs must be long or weakly informative. On the flu graph (below) the dates of the two variants differ by less than 0.09 years at 99% of nodes
- **Hidden nodes are not dated by the data**: their coverage is 0.5 to 0.6. Optimize holds the edge below such a node at length 0 (`TimelessNode`, [kb/algo/optimization.md](../algo/optimization.md)), so the time passes place it at its child, with a region as narrow as the child's. In the full trees they also pull the neighbouring deep nodes (RMSE 0.6 to 1.3 years against 0.2 to 0.64 on the collapsed trees, from a few nodes near the roots). Real inputs come from tree inference and do not contain them
- **Refinement rounds are needed**: optimize fixes only the sum of the two edges below a segment root. With a single time pass (before commit `09294169f`), the arbitrary split moved deep nodes of the full trees by up to 21 years (high, rep02). With clock branch lengths fed back into the marginal reconstruction, the largest error at a visible node is 12.3 years on the full trees and 6.3 years on the collapsed trees (high), 8.0 and 1.5 years (low); the large errors sit at the oldest nodes, near the roots of the graph

## Flu (four segments, 2839 nodes)

`timetree --segment-trees` on the named TreeKnit trees of [kb/reports/2026-11-15-treeknit-output-to-named-segment-trees.md](2026-11-15-treeknit-output-to-named-segment-trees.md) at `--clock-rate 0.003` with `--confidence`: 43 s and 350 MB with cross-parent messages (optimize about 15 s; 3 time rounds with 9, 20 and 16 sweeps; the second reached the sweep limit at a relative change of 2 x 10^-4), 19 s and 310 MB with `--no-cross-parent-messages`.

## Open

- **Default of the cross-parent messages**: exact tests favour them; on simulated and flu graphs they change nothing measurable and double the run time
- **Hidden nodes**: contract them for the time inference, or report them without a date
- **Sweep tolerance**: the flu graph needs more than 20 sweeps in one round to reach 10^-4 of the date span
- Clock rate estimation, coalescent prior, rooting and polytomy resolution on segment graphs: [kb/issues/M-timetree-segment-graphs-fixed-clock-only.md](../issues/M-timetree-segment-graphs-fixed-clock-only.md)
