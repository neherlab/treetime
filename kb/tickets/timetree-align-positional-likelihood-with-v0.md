# Align positional likelihood metric with v0

v1 uses a different formula for the positional log-likelihood convergence metric than v0. Both metrics track the same quantity (how well inferred times explain branch-length distributions) and trend in the same direction during convergence, but produce different numerical values.

- v1: `compute_positional_log_lh()` (`#compute_positional_log_lh`) at
  [`packages/treetime/src/timetree/convergence/likelihood.rs#L35-L76`](../../packages/treetime/src/timetree/convergence/likelihood.rs#L35-L76)
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

`log_lh_pos` is a reported diagnostic only: the convergence gate `has_converged()` depends solely on `max_time_change` and `n_resolved` (`packages/treetime/src/timetree/convergence/metrics.rs:49-54`), not on any log-likelihood. Correcting the metric changes reported `log_lh_pos`/`log_lh_total`, not inferred dates or the stopping iteration.

## Correctness fix required independently of parity

The v1 per-edge term is also internally wrong, separate from the v0-parity question. `Distribution::<NegLog>::eval` returns the neg-log ordinate `y = -ln(p / p_peak) >= 0` (`packages/treetime-distribution/src/distribution_core/distribution.rs:146`), not a probability. v1 sums `y.ln()`; the documented log-probability is `-y` (v0 sums `-BLI = -y`). Fix the accumulation to a genuine log-probability and rework the `p > 0.0` guard and "zero or negative probability" message (`packages/treetime/src/timetree/convergence/likelihood.rs:56-67`), which assume a plain probability and wrongly drop the mode (`y = 0`). This fix stands regardless of whether the calendar-time-vs-substitution and root-term differences are aligned to v0.

## Related issues

- Source: [kb/issues/M-timetree-positional-likelihood-metric.md](../issues/M-timetree-positional-likelihood-metric.md) -- delete after full resolution
