# Timetree on segment graphs uses loopy message passing, with the subtree alone sent to each parent by default

A reassortant node of a segment graph has several parents, so the time graph has undirected cycles, and the two passes of the tree no longer give exact marginals. v1 keeps the passes node-local, as loopy belief propagation (<a id="cite-1"></a>[Koller and Friedman 2009](https://mitpress.mit.edu/9780262013192/probabilistic-graphical-models/) [[1](#ref-1)]):

- **Backward**: a node with several parents sends each parent its subtree distribution convolved across the edge. With `--cross-parent-messages` it multiplies in the downward messages that the other parents sent in the previous sweep
- **Forward**: a node multiplies its subtree distribution with the downward message of each parent, and its committed time is at least the latest committed time of its parents
- **Sweeps**: backward and forward passes repeat until the largest change of a likely node time is at most 10^-4 of the date span, at most 20 sweeps. Without cross-parent messages one sweep is final, because the backward pass then reads nothing from the previous sweep

All segments share one clock rate, given with `--clock-rate` ([kb/issues/M-timetree-segment-graphs-fixed-clock-only.md](../issues/M-timetree-segment-graphs-fixed-clock-only.md)).

**Type**: Approximation of the time posterior on graphs with cycles (new v1 model). Approved.

Code: `fn propagate_distributions_backward_on_graph` [packages/treetime/src/timetree/inference/backward_pass.rs](../../packages/treetime/src/timetree/inference/backward_pass.rs), `fn propagate_distributions_forward_on_graph` [packages/treetime/src/timetree/inference/forward_pass.rs](../../packages/treetime/src/timetree/inference/forward_pass.rs), `fn sweep_time_passes` [packages/treetime/src/timetree/inference/runner.rs](../../packages/treetime/src/timetree/inference/runner.rs). Algorithm: [kb/algo/timetree.md](../algo/timetree.md#graphs-with-several-parents).

## Accuracy

Maximum CDF error against the exact posterior on a grid ([`packages/treetime/src/timetree/inference/__tests__/test_reassortment_passes.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_reassortment_passes.rs), values from [kb/algo/timetree.md](../algo/timetree.md#graphs-with-several-parents)):

| graph                  | with cross-parent messages | subtree alone (default) |
| ---------------------- | -------------------------- | ----------------------- |
| polytree (no cycle)    | about 10^-4, exact after two sweeps | more than 4 times larger at a parent |
| one cycle              | 0.003 to 0.005             | 0.05 to 0.08            |

On simulated segment graphs and on the H3N2 graph the two variants agree ([kb/reports/2026-11-15-segment-validation-timetree.md](../reports/2026-11-15-segment-validation-timetree.md)):

- **Simulated collapsed segment trees**: RMSE within 0.01 years of each other; the 90% regions cover 85 to 91% of the true dates, at nodes with one parent and with several parents alike
- **H3N2, 4 segments, 2839 nodes**: the dates differ by less than 0.09 years at 99% of nodes; 19 s without cross-parent messages, 43 s with them, where one refinement round reached the limit of 20 sweeps

## Why the subtree alone is the default

- **No measured gain**: the cross-parent messages change no date by a measurable amount on the simulated graphs or on H3N2
- **Cost**: they need up to 20 sweeps per refinement round in place of one
- **Known error**: without them, an ancestor above a cycle receives the evidence of the subtree below the cycle once along each parent path, which gives the larger CDF error on small cycles

The log marginal likelihood that the passes report sums over the roots of the graph and is not exact on a graph with cycles, with or without cross-parent messages.

## Options not chosen

- **Conditioning on one parent edge per reassortant node** and combining the runs, as v0 does with its focal trees: the runs disagree on shared nodes ([segment-graph-replaces-v0-arg-masks.md](segment-graph-replaces-v0-arg-masks.md))
- **Exact treatment of the cycles**, for example a junction tree over the nodes of each cycle: its cost grows with the size of the cycles

## References

1. <a id="ref-1"></a> Koller, Daphne, and Nir Friedman. 2009. _Probabilistic Graphical Models: Principles and Techniques._ MIT Press. ISBN 978-0-262-01319-2. https://mitpress.mit.edu/9780262013192/probabilistic-graphical-models/ [↩](#cite-1)
