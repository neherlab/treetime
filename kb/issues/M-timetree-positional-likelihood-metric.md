# Positional likelihood metric differs from v0

v1 uses a different formula for the positional likelihood convergence metric than v0. Both metrics track the same quantity (how well inferred times explain branch-length distributions) and trend in the same direction during convergence, but produce different numerical values.

- v1: `compute_positional_likelihood()` (`#compute_positional_likelihood`) at
  [`packages/treetime/src/commands/timetree/convergence/likelihood.rs#L35-L76`](../../packages/treetime/src/commands/timetree/convergence/likelihood.rs#L35-L76)
- v0: `timetree_likelihood()` (`#timetree_likelihood`) at
  [`packages/legacy/treetime/treetime/clock_tree.py#L621-L639`](../../packages/legacy/treetime/treetime/clock_tree.py#L621-L639)
- v0 usage: tracked as `positional_LH` in the convergence log at
  [`packages/legacy/treetime/treetime/treetime.py#L274`](../../packages/legacy/treetime/treetime/treetime.py#L274) and
  [`packages/legacy/treetime/treetime/treetime.py#L357`](../../packages/legacy/treetime/treetime/treetime.py#L357)

## v0 computation

v0 has two paths depending on `time_marginal`:

- Marginal: integrates the root's `marginal_pos_LH` distribution over its
  support (a single scalar from numerical integration of the root node's time-position distribution).
- Joint (default): sums `-BLI(branch_length)` over all non-root edges, where
  `BLI` is the branch length interpolator (returns neg-log probability evaluated at evolutionary distance in substitutions/site). Adds the root sequence log-likelihood.

## v1 computation

v1 sums `ln(dist.eval(child_time - parent_time))` over all edges with assigned times, where `dist` is the branch-length distribution evaluated at the calendar time duration (not evolutionary distance). No root sequence term.

## Differences

- v0 evaluates at evolutionary distance (`branch_length` in subs/site); v1
  evaluates at calendar time duration (`child_time - parent_time`).
- v0 joint path includes a root sequence log-likelihood term; v1 does not.
- v0 marginal path integrates the root distribution; v1 has no marginal path
  for this metric.

The metric drives convergence detection in the EM refinement loop. Different numerical values can cause the loop to terminate at a different iteration, producing different final timetree results.

This is a candidate for alignment with v0's formula if parity is desired, or for explicit documentation as a v1-specific improvement if the new metric is scientifically preferable. Needs discussion.
