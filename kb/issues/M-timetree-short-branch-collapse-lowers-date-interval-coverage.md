# The pre-loop short-branch collapse lowers the coverage of the date intervals

> [!IMPORTANT]
> **Investigation required.** Find a collapse that keeps the date intervals calibrated, and the cause of the lost coverage. Evidence to collect: the coverage of the date intervals of the surviving nodes against the coverage of the same nodes without the collapse, per node depth and per number of children that the node gains, and the same measurement on v0.

## Problem

Before the timetree loop, `timetree` collapses each non-root internal edge shorter than $0.1/L$ whose states equal the states of its parent, as v0 does (`fn collapse_short_internal_edges()` in [packages/treetime/src/optimize/topology/collapse.rs](../../packages/treetime/src/optimize/topology/collapse.rs), called from [packages/treetime/src/timetree/pre_loop.rs](../../packages/treetime/src/timetree/pre_loop.rs)). The collapse is on by default; `--no-collapse-short-branches` turns it off.

On simulated data the collapse makes the 90% date intervals of the internal nodes too narrow. The test `test_sim_timetree_date_intervals_cover_the_true_dates` ([packages/treetime/src/timetree/__tests__/test_sim_timetree_recovery.rs](../../packages/treetime/src/timetree/__tests__/test_sim_timetree_recovery.rs)) runs the standard simulated cell at the true $T_c$ and clock rate, and requires the fraction of true internal node dates inside their interval to be $0.9 \pm 0.03$:

- With the collapse, paired by clade over the internal nodes that survive it: coverage 0.834; error mean -0.065527, standard deviation 0.114811, n 20446, interval of the mean [-0.068008, -0.063045]. The test fails
- Without the collapse: the test passes
- In the standard cell the collapse removes 8 of the 19 internal nodes

The validation suite shows the same loss ([kb/reports/2026-11-08-timetree-short-branch-collapse-on-true-and-reconstructed-topologies.md](../reports/2026-11-08-timetree-short-branch-collapse-on-true-and-reconstructed-topologies.md)). Skyline prior with `--confidence`, estimated clock rate, 400 tips, collapse / no collapse:

- True topology, every surviving internal node: 0.745 / 0.900 at $L\mu=1$, 0.793 / 0.903 at $L\mu=3$. The intervals with the collapse are wider on average (0.50 / 0.35 y at $L\mu=1$), so the dates of the polytomies are biased, not only too narrow
- FastTree topology, inferred nodes that are true clades: 0.73 / 0.80 at $L\mu=1$, 0.76 / 0.87 at 3, 0.81 / 0.87 at 10. At $L\mu=10$: 0.78 / 0.88 at 1200 tips, 0.75 / 0.87 at 4000 tips

On the FastTree topology the collapse still gives the more accurate dates at $L\mu \le 3$ (400 tips) and up to $L\mu=30$ (4000 tips), so turning it off is no fix. A candidate is to widen the interval of a polytomy, see the heuristic in [kb/issues/M-timetree-date-intervals-ignore-topology-uncertainty.md](M-timetree-date-intervals-ignore-topology-uncertainty.md#fix-direction).

The simulation tests therefore date the true topology: the simulated cell sets `TimetreeParams::collapse_short_branches` to `false` (`fn timetree_cell_params()` in [packages/treetime/src/test_utils/simulation.rs](../../packages/treetime/src/test_utils/simulation.rs)), so they do not cover the calibration of a collapsed tree.

## Fix direction

Open: a collapse that keeps the date intervals calibrated.

## Validation

`test_sim_timetree_date_intervals_cover_the_true_dates` passes with the collapse on.

## Proposed implementation plan

The cited locations are current (`fn collapse_short_internal_edges()`, its call in `pre_loop.rs`, `fn timetree_cell_params()` at [packages/treetime/src/test_utils/simulation.rs#L67](../../packages/treetime/src/test_utils/simulation.rs#L67) with `collapse_short_branches: false`, and the `--no-collapse-short-branches` flag in `packages/app-commands/src/commands/timetree/args.rs`).

1. Collect the evidence of the Investigation required block: in `test_sim_timetree_date_intervals_cover_the_true_dates`, add a measurement mode with the collapse on that splits the coverage and the signed error of the surviving internal nodes by node depth and by the number of children the node gains. Run the same simulated cells through v0 with the true topology and its own collapse, and compare
2. Locate the cause from the split: when the loss concentrates at nodes that gain children, the interval of a polytomy misses the time spread of the merged nodes; when it concentrates at their parents, the collapse moves the branch-length evidence of the removed edges
3. Implement the fix at the located cause, in the owner of the interval (`packages/treetime/src/timetree/confidence.rs`) or of the collapse (`packages/treetime/src/optimize/topology/collapse.rs`). For the polytomy case, the candidate is the interval widening of [kb/issues/M-timetree-date-intervals-ignore-topology-uncertainty.md](M-timetree-date-intervals-ignore-topology-uncertainty.md#fix-direction); a scientific change of the intervals needs a user-approved entry in `kb/decisions/`
4. Set `collapse_short_branches: true` in `fn timetree_cell_params()`, so the simulation tests cover the collapsed tree
5. Validate: `just test-slow-rs` (on user request) with `test_sim_timetree_date_intervals_cover_the_true_dates` at $0.9 \pm 0.03$, the oracle being the simulated true dates, and the validation suite of the report on true and FastTree topologies, where the dates with the collapse must stay at least as accurate as without it
