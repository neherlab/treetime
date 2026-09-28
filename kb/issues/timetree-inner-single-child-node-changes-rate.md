# An undated single-child node inside the tree changes the timetree result

## Problem

Adding an undated single-child node on a branch changes the clock rate that `treetime` (timetree) reports. In theory the node carries no data, so the result should not change: substitution probabilities satisfy $P(t_1)P(t_2) = P(t_1 + t_2)$, and root-to-tip distances and variances add along a path.

`resolve_polytomies()` removes undated single-child nodes, but only after the first timetree iterations.

## Evidence

Inputs and commands: runs `no-stem` and `inner-node` in [kb/reports/reroot-single-child-root.md](../reports/reroot-single-child-root.md#reproduction). The two trees differ only by the undated node `u1` above `i359805`.

- Rate `8.452e-04` without `u1` and `8.457e-04` with it
- The difference appears in joint (`--time-marginal never`) and marginal (`--time-marginal always`) runs
- `treetime clock` gives the same rate and root for both trees, so the difference comes from the timetree iterations, not from the initial regression and reroot

## Open question

The cause is not traced. Candidates to check: the per-branch `clock_length` used in the covariance of the regression after the first iteration, and the discretized convolution of branch-length distributions across the extra node.
