# Relaxed-clock rates are reset in every round that resolves a polytomy

In each refinement round, v1 fits the relaxed-clock rate multipliers (gamma) and then resolves polytomies. When polytomy resolution changes the topology, v1 resets gamma to 1.0 on every edge, so that round's time inference, clock branch lengths and clock update run without the relaxed clock.

## v1 behavior

- `Refinement::run` calls `apply_relaxed_clock` before `refine_topology` ([packages/treetime/src/timetree/refinement.rs](../../packages/treetime/src/timetree/refinement.rs))
- `refine_topology` calls `TimetreeState::reset_date_edges_for_topology_change`, which sets `gamma = 1.0` on every edge ([packages/treetime/src/timetree/timetree_state.rs](../../packages/treetime/src/timetree/timetree_state.rs))
- New edges already default to 1.0 without the reset, and re-parented children keep their edge keys

## v0 behavior

- Same order within a round: relaxed clock, then polytomies ([packages/legacy/treetime/treetime/treetime.py#L311-L343](../../packages/legacy/treetime/treetime/treetime.py#L311-L343))
- After resolution, `make_time_tree` rebuilds the branch-length interpolators but copies each existing node's gamma first ([packages/legacy/treetime/treetime/clock_tree.py#L349-L370](../../packages/legacy/treetime/treetime/clock_tree.py#L349-L370)). Only new nodes start at 1.0 ([packages/legacy/treetime/treetime/branch_len_interpolator.py#L29](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L29))

## Impact

With `--relax` and `--resolve-polytomies`, the relaxed clock has no effect in rounds that resolve polytomies. v0 parity is: keep gamma on existing edges, 1.0 on new edges.
