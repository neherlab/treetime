# Mass-sized time distributions produce node times that break downstream invariants

## Symptom and reproduction

With probability-mass grid sizing enabled for timetree distributions (`rewindow_to_mass` at branch-length construction and after the backward/forward composite steps) and the forward-pass division bounded to the divisor's real support (so it no longer collapses to `Empty`), the forward pass completes but assigns node times that violate two downstream invariants on some datasets:

- **Inverted branch.** A committed child time can precede its parent's, so the branch duration is negative and lies outside its branch-length distribution's support (past the hard left boundary at `t = 0`).
- **Seeded polytomy resolution flips.** `resolve_polytomies` reaches a different topology outcome than before mass sizing: a polytomy that used to be resolved is now left unresolved. The resolver (`simulate_subtree`) is seeded-stochastic and sensitive to the exact child/parent times, so a node-time shift changes the sampled history and the resolution decision.

`test_round_rebuilds_complete_coalescent_state_after_topology_change` and `test_round_unchanged_topology_recomputes_missing_time` (`packages/treetime/src/timetree/__tests__/test_round.rs`) fail and are `#[ignore]`d: polytomy resolution merges no lineage (`Changed { resolved_nodes: 1 }` -> `Unchanged`), which cascades to `Polytomy resolution requires an inferred time for node N, but it has none` on the follow-up refinement.

The three pipeline tests in `packages/app-commands/src/commands/timetree/__tests__/test_pipeline.rs` that this issue also blocked failed only on the positional log-likelihood metric; the metric was removed in favour of `log_lh_time` and they run again.

## Impact and scope

- The inverted branch is a genuine correctness symptom: a committed node time is inconsistent with its incident branch-length support.
- The polytomy `resolved_nodes` assertions predate mass sizing (the value was asserted and passing before the mass-sizing collapse ignore was added), so they encode pre-mass-sizing node times.

## Root cause

Mass sizing stores each distribution on its own mass-bounded domain instead of the heuristic peak-multiple grid, which moves the stored grid extent and can shift `likely_time` (the peak) by up to a grid spacing. Downstream consumers assume the pre-mass-sizing times; `resolve_polytomies` is a seeded sampler whose output is not stable under sub-grid time perturbations.

## Confirmed mechanism (topology violation, not clipping)

On `test_pipeline_timetree_convergence` the offending edge was `17 -> 9`: `parent_time = 2003.00717`, `child_time = 2003.00274`, so `child_time - parent_time = -4.43e-3` (and `-7.64e-3` on the next iteration), against a branch-length support of `(0.0, 5.87)`.

This is mechanism (b): a marginal node-time inversion (child before parent), overlapping [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md). It is **not** mechanism (a) clipping: the duration is not a clipped positive tail, it is negative. Mass sizing re-grids the parent posterior and shifts its committed `likely_time` by a sub-grid amount (~1.6 days here), which is enough to push an already near-zero branch past its earlier, near-fixed child. `committed_time` projects an inferred node up to its parent but does not prevent an inferred parent from landing after an earlier child, so the inversion survives.

## Fix approach

- Resolve the marginal node-time contract in [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md). This decides whether an inferred parent may land after an earlier child and is the substantive fix for the inverted branch. Decision-required.
- The seeded polytomy flip is a separate consequence of the same node-time perturbation: `resolve_polytomies` (`simulate_subtree`) is sensitive to sub-grid time shifts. Making mass sizing mode-preserving (committed `likely_time` matches the pre-mass-sizing peak within tolerance) would keep the resolver stable; the resolver's sensitivity to sub-grid perturbation is otherwise its own robustness concern.

Any fix must keep the golden-master reference runs unchanged and restore the two ignored round tests to passing.

## Related

- [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md): marginal node times can invert parent-child order. This is the confirmed root cause of the inverted branch above; its topology contract decision is the substantive fix.
- [M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md): branch-grid resolution concern for the same distributions.
