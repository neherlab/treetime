# Positional likelihood metric differs from v0

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

## v1 formula is internally incorrect, not just a different valid metric

Beyond the v0/v1 divergence above, the v1 per-edge term is mathematically wrong. `Distribution::<NegLog>::eval` returns the stored neg-log ordinate `y = -ln(p / p_peak) >= 0` (`packages/treetime-distribution/src/distribution_core/distribution.rs:146`, `eval = f.interp`), not a plain probability. v1 sums `y.ln()`; the log-probability it is documented to sum is `-y`.

- v0's joint path sums `-BLI(branch_length)`, i.e. `-y`, the log-probability. v1 sums `+ln(y)`, an unrelated quantity.
- The direction is inverted across edges: as an edge fit worsens, `y` grows, the true log-probability `-y` decreases, but the summed `ln(y)` increases.
- The `p > 0.0` guard and the "zero or negative probability" debug message (`packages/treetime/src/timetree/convergence/likelihood.rs:56-67`) assume a plain probability. The actual mode (`y = 0`, highest probability) fails `p > 0.0`, is logged as "zero or negative probability", and is dropped from the sum.

So this is not a "possibly preferable v1 metric pending a parity discussion": the v1 formula must be corrected to a genuine log-probability (`-y`, or evaluate in plain space and take its log) regardless of whether the calendar-time-vs-substitution and root-term differences above are aligned to v0.

## Impact

`log_lh_pos` feeds `log_lh_total` and the `ConvergenceMetrics` written to the tracelog and the per-iteration info log (`packages/treetime/src/timetree/convergence/optimizer.rs:59-91`). It is a reported diagnostic only: `has_converged()` depends solely on `max_time_change` and `n_resolved` (`packages/treetime/src/timetree/convergence/metrics.rs:49-54`), so the defect does not change inferred node dates, tree topology, or the convergence stopping decision. The observable effect is an incorrect `log_lh_pos` (and therefore `log_lh_total`) surfaced to the user.

(An earlier statement that this metric "drives convergence detection" does not hold for the current convergence gate.)

## Interaction with the branch-length right-boundary tail

The branch-length distribution now carries a soft `Linear` right tail (`packages/treetime/src/timetree/inference/branch_length_likelihood.rs`), so an edge whose inferred duration exceeds the grid's `t_max` is evaluable instead of failing `eval`. Such edges previously hit the `Err` arm and were dropped from this sum; they now flow through the incorrect `ln(y)` accumulation above, so more edges are affected than before.

This is a candidate for alignment with v0's formula if parity is desired, but the internal-correctness fix (`ln(y)` to `-y`) is required independently of that decision.

## Related tickets

- [kb/tickets/timetree-align-positional-likelihood-with-v0.md](../tickets/timetree-align-positional-likelihood-with-v0.md)
