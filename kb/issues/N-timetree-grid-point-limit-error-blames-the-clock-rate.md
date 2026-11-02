# The grid point limit error of timetree blames the clock rate whatever the cause

When the time inference requests a grid above `--max-grid-points`, `fn explain_grid_point_limit()` [packages/treetime/src/timetree/inference/runner.rs#L155-L184](../../packages/treetime/src/timetree/inference/runner.rs#L155-L184) wraps the error in a cause that names the clock rate for every run:

- **Fixed rate**: "The fixed --clock-rate R may be far from the rate the data supports, which makes branch-time distributions very narrow. Check --clock-rate"
- **Estimated rate**: "The estimated clock rate R may be far from the true rate, which makes branch-time distributions very narrow; errors in the input dates and weak temporal signal cause such estimates. Check the input dates, set --clock-rate"

The message leads users to change a correct rate. In the known cases, the rate is often right and another mechanism requests the grid:

- **Chains of zero-length branches**: rubella E1, oropouche S, and the two H5N1 cattle-outbreak builds of Nextstrain stop with the fixed-rate message at the clock rates of their workflows. Each finishes at the same rate once its zero-length branches are collapsed ([H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md](H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md))
- **Near-flat soft tail**: the Nextstrain OC43 build stops with the fixed-rate message on a branch grid over 26 180 years, and finishes at the same rate with `--no-indels` ([M-timetree-near-flat-soft-tail-stretches-the-time-grid.md](M-timetree-near-flat-soft-tail-stretches-the-time-grid.md))
- **Wrong coalescent time scale**: on the standard simulated cell (`fn timetree_cell()` [packages/treetime/src/test_utils/simulation.rs#L49](../../packages/treetime/src/test_utils/simulation.rs#L49), true clock rate fixed), runs with a wrong fixed coalescent time scale $T_c$ such as 0.5 or 1.0, against a true $T_c$ of 5, sometimes stop with the fixed-rate message, although the fixed rate is the true one

The message fits the case it was written for: on `data/lassa/L/20`, with no temporal signal ($r = 0.023$), the estimated rate is the cause ([M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md](M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md)).

## Fix direction

Report what the run knows about the request: the operation and the node or edge that requests the grid, its span and spacing, and which operand sets each. Name the clock rate only when the narrow operand is a branch distribution whose width the rate sets.

## Validation

- A test that exceeds the limit through a narrow node posterior gets a message that names the node and does not name the clock rate
- The `grid-point-limit` variant of `dev/smoke.toml` (mpox at `--clock-rate=0.5`) still names the clock rate
