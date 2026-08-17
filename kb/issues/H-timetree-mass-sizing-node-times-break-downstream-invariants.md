# Mass-sized time distributions produce node times that break downstream invariants

## Symptom and reproduction

With probability-mass grid sizing enabled for timetree distributions (`rewindow_to_mass` at branch-length construction and after the backward/forward composite steps) and the forward-pass division bounded to the divisor's real support (so it no longer collapses to `Empty`), the forward pass completes but assigns node times that violate two downstream invariants on some datasets:

- **Positional log-likelihood `+inf`.** `compute_positional_log_lh` evaluates each edge's branch-length distribution at `child_time - parent_time` and reads zero probability there (`-ln 0 = +inf`), so the pipeline's `log_lh_pos < 0` invariant fails. The inferred branch duration lands where the branch-length distribution carries no probability (the `HardApproach` left tail past `t_hard`, i.e. a near-zero or negative duration).
- **Seeded polytomy resolution flips.** `resolve_polytomies` reaches a different topology outcome than before mass sizing: a polytomy that used to be resolved is now left unresolved. The resolver (`simulate_subtree`) is seeded-stochastic and sensitive to the exact child/parent times, so a node-time shift changes the sampled history and the resolution decision.

Reproduce on `dev` with the division bound in place:

- `test_pipeline_timetree_convergence` (`packages/treetime/src/commands/timetree/__tests__/test_pipeline.rs`) fails with `Positional log-likelihood must be negative, got inf`.
- `test_refinement_rebuilds_complete_coalescent_state_after_topology_change` and `test_refinement_unchanged_topology_recomputes_missing_time` (`packages/treetime/src/timetree/__tests__/test_refinement.rs`) fail: the polytomy resolves 0 nodes (`Changed { resolved_nodes: 1 }` -> `Unchanged`), which cascades to `Polytomy resolution requires an inferred time for node N, but it has none` on the follow-up refinement.

## Impact and scope

- Blocks restoring the four integration tests that the forward-pass division collapse originally forced `#[ignore]` (see [H-timetree-mass-sizing-collapses-forward-pass-division.md](H-timetree-mass-sizing-collapses-forward-pass-division.md)). Bounding the division removes the collapse but exposes this next layer.
- The positional log-likelihood `+inf` is a genuine correctness symptom: a committed node time is inconsistent with its incident branch-length support.
- The polytomy `resolved_nodes` assertions predate mass sizing (the value was asserted and passing before the mass-sizing collapse ignore was added), so they encode pre-mass-sizing node times.

## Root cause (partial)

Mass sizing stores each distribution on its own mass-bounded domain instead of the heuristic peak-multiple grid, which moves the stored grid extent and can shift `likely_time` (the peak) by up to a grid spacing. Downstream consumers assume the pre-mass-sizing times:

- `compute_positional_log_lh` assumes every committed duration lands inside its branch-length distribution's support.
- `resolve_polytomies` is a seeded sampler whose output is not stable under sub-grid time perturbations.

Two candidate mechanisms for the positional `+inf` are not yet distinguished:

- A near-zero or negative duration (a node timed at or before its parent), overlapping [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md).
- A legitimate positive duration clipped out of a too-narrow branch-length mass domain (the construction pilot grid clipping a real tail before the mass rewindow measures it).

## Fix approach

Not decided. Options:

- Determine which mechanism produces the positional `+inf` (instrument the edge and duration that reads zero probability). If it is a clipped branch-length mass domain, make the construction pilot grid mass-sufficient so the rewindow never measures a clipped tail. If it is a topology violation, resolve the marginal node-time monotonicity separately.
- Reconsider whether message-passing posteriors should be stored mass-sized at all, versus mass sizing only final per-node outputs in a way that never feeds a subsequent division or breaks branch-support consistency.
- Make mass sizing mode-preserving so the committed `likely_time` matches the pre-mass-sizing peak within tolerance, keeping the seeded polytomy resolution stable.

Any fix must keep the golden-master reference runs unchanged and restore the four ignored integration tests to passing.

## Related

- [H-timetree-mass-sizing-collapses-forward-pass-division.md](H-timetree-mass-sizing-collapses-forward-pass-division.md): the division-sampling cause, now bounded so the division no longer collapses.
- [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md): marginal node times can invert parent-child order.
- [M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md): branch-grid resolution concern for the same distributions.
