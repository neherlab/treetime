# The pre-loop short-branch collapse lowers the coverage of the date intervals

> [!IMPORTANT]
> **Investigation required.** Find a collapse that keeps the date intervals calibrated, and the cause of the lost coverage. Evidence to collect: the coverage of the date intervals of the surviving nodes against the coverage of the same nodes without the collapse, per node depth and per number of children that the node gains, and the same measurement on v0.

## Problem

Before the timetree loop, `timetree` collapses each non-root internal edge shorter than $0.1/L$ whose states equal the states of its parent, as v0 does (`fn collapse_short_internal_edges()` in [packages/treetime/src/optimize/topology/collapse.rs](../../packages/treetime/src/optimize/topology/collapse.rs), called from [packages/treetime/src/timetree/pre_loop.rs](../../packages/treetime/src/timetree/pre_loop.rs)). The collapse is on in every CLI run, and no flag turns it off.

On simulated data the collapse makes the 90% date intervals of the internal nodes too narrow. The test `test_sim_timetree_date_intervals_cover_the_true_dates` ([packages/treetime/src/timetree/__tests__/test_sim_timetree_recovery.rs](../../packages/treetime/src/timetree/__tests__/test_sim_timetree_recovery.rs)) runs the standard simulated cell at the true $T_c$ and clock rate, and requires the fraction of true internal node dates inside their interval to be $0.9 \pm 0.03$:

- With the collapse, paired by clade over the internal nodes that survive it: coverage 0.834; error mean -0.065527, standard deviation 0.114811, n 20446, interval of the mean [-0.068008, -0.063045]. The test fails
- Without the collapse: the test passes
- In the standard cell the collapse removes 8 of the 19 internal nodes

The simulation tests therefore date the true topology: the simulated cell sets `TimetreeParams::collapse_short_branches` to `false` (`fn timetree_cell_params()` in [packages/treetime/src/test_utils/simulation.rs](../../packages/treetime/src/test_utils/simulation.rs)), so they do not cover the calibration of a collapsed tree.

## Fix direction

Open: a collapse that keeps the date intervals calibrated.

## Validation

`test_sim_timetree_date_intervals_cover_the_true_dates` passes with the collapse on.
