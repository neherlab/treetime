# Mass-sized time distributions produce node times that break downstream invariants

## Symptom and reproduction

With probability-mass grid sizing enabled for timetree distributions (`rewindow_to_mass` at branch-length construction and after the backward/forward composite steps) and the forward-pass division bounded to the divisor's real support (so it no longer collapses to `Empty`), the forward pass completes but assigns node times that violate two downstream invariants on some datasets:

- **Positional log-likelihood `+inf`.** `compute_positional_log_lh` evaluates each edge's branch-length distribution at `child_time - parent_time` and reads zero probability there (`-ln 0 = +inf`), so the pipeline's `log_lh_pos < 0` invariant fails. The inferred branch duration lands where the branch-length distribution carries no probability (the `HardApproach` left tail past `t_hard`, i.e. a near-zero or negative duration).
- **Seeded polytomy resolution flips.** `resolve_polytomies` reaches a different topology outcome than before mass sizing: a polytomy that used to be resolved is now left unresolved. The resolver (`simulate_subtree`) is seeded-stochastic and sensitive to the exact child/parent times, so a node-time shift changes the sampled history and the resolution decision.

Reproduce on `dev` with the division bound in place:

- `test_pipeline_timetree_convergence` (`packages/treetime/src/commands/timetree/__tests__/test_pipeline.rs`) fails with `Positional log-likelihood must be negative, got inf`.
- `test_refinement_rebuilds_complete_coalescent_state_after_topology_change` and `test_refinement_unchanged_topology_recomputes_missing_time` (`packages/treetime/src/timetree/__tests__/test_refinement.rs`) fail: the polytomy resolves 0 nodes (`Changed { resolved_nodes: 1 }` -> `Unchanged`), which cascades to `Polytomy resolution requires an inferred time for node N, but it has none` on the follow-up refinement.

## Impact and scope

- Blocks restoring the four integration tests that the now-fixed forward-pass division collapse originally forced `#[ignore]`. Bounding the division quotient to the divisor's grid support (`divisor_tail_extends_quotient` in `packages/treetime-distribution/src/distribution_ops/divide.rs`) removed the empty-distribution collapse but exposed this next layer.
- The positional log-likelihood `+inf` is a genuine correctness symptom: a committed node time is inconsistent with its incident branch-length support.
- The polytomy `resolved_nodes` assertions predate mass sizing (the value was asserted and passing before the mass-sizing collapse ignore was added), so they encode pre-mass-sizing node times.

## Root cause (partial)

Mass sizing stores each distribution on its own mass-bounded domain instead of the heuristic peak-multiple grid, which moves the stored grid extent and can shift `likely_time` (the peak) by up to a grid spacing. Downstream consumers assume the pre-mass-sizing times:

- `compute_positional_log_lh` assumes every committed duration lands inside its branch-length distribution's support.
- `resolve_polytomies` is a seeded sampler whose output is not stable under sub-grid time perturbations.

## Confirmed mechanism (topology violation, not clipping)

Instrumenting `compute_positional_log_lh` on `test_pipeline_timetree_convergence` identified the offending edge as `17 -> 9`: `parent_time = 2003.00717`, `child_time = 2003.00274`, so `time_diff = child_time - parent_time = -4.43e-3` (and `-7.64e-3` on the next iteration). The duration is negative (the child is timed before its parent), and the branch-length distribution support is `(0.0, 5.87)`, so the negative duration falls beyond the hard left boundary at `t = 0`, where `interp` returns `Y::probability_zero()` = `+inf` under `NegLog`.

This is mechanism (b): a marginal node-time inversion (child before parent), overlapping [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md). It is **not** mechanism (a) clipping: the duration is not a clipped positive tail, it is genuinely negative. Mass sizing re-grids the parent posterior and shifts its committed `likely_time` by a sub-grid amount (~1.6 days here), which is enough to push an already near-zero branch past its earlier, near-fixed child. `set_likely_time` projects an inferred node up to its parent but does not prevent an inferred parent from landing after an earlier child, so the inversion survives.

The `+inf` panic (rather than a benign `-inf`) is compounded by a unit bug in the metric: `compute_positional_log_lh` (`packages/treetime/src/timetree/convergence/likelihood.rs`) evaluates the `NegLog` branch-length distribution and treats the returned neg-log ordinate as a plain probability (`p > 0.0` then `total += p.ln()`). The log-probability contribution should be `-p` (`= log prob`). For a valid duration the wrong formula happens to stay negative; for the beyond-boundary `+inf` neg-log it yields `+inf` instead of the semantically correct `-inf`, tripping the `log_lh_pos < 0` assertion with a misleading sign. This metric defect is independent of the node-time inversion and is tracked separately (see Related).

## Fix approach

Not decided. Options:

- Determine which mechanism produces the positional `+inf` (instrument the edge and duration that reads zero probability). If it is a clipped branch-length mass domain, make the construction pilot grid mass-sufficient so the rewindow never measures a clipped tail. If it is a topology violation, resolve the marginal node-time monotonicity separately.
- Reconsider whether message-passing posteriors should be stored mass-sized at all, versus mass sizing only final per-node outputs in a way that never feeds a subsequent division or breaks branch-support consistency.
- Make mass sizing mode-preserving so the committed `likely_time` matches the pre-mass-sizing peak within tolerance, keeping the seeded polytomy resolution stable.

Any fix must keep the golden-master reference runs unchanged and restore the four ignored integration tests to passing.

## Related

- [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md): marginal node times can invert parent-child order. This is the confirmed root cause of the positional `+inf` above; its topology contract decision is the substantive fix.
- [M-timetree-positional-likelihood-metric.md](M-timetree-positional-likelihood-metric.md): the `eval` returns a neg-log ordinate but the metric sums `y.ln()` instead of `-y`. This unit bug turns the beyond-boundary `+inf` neg-log into a `+inf` sum (crash) instead of the correct `-inf`.
- [M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md): branch-grid resolution concern for the same distributions.
