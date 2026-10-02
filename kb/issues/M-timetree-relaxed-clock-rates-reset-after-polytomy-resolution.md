# Relaxed-clock rates are reset in every round that resolves a polytomy

In each refinement round, v1 fits the relaxed-clock rate multipliers (gamma) and then resolves polytomies. When polytomy resolution changes the topology, v1 resets gamma to 1.0 on every edge, so that round's time inference, clock branch lengths and clock update run without the relaxed clock.

## v1 behavior

- `refinement_round` calls `relax_clock` (which calls `apply_relaxed_clock`) before `refine_topology` ([packages/treetime/src/timetree/round.rs](../../packages/treetime/src/timetree/round.rs))
- `refine_topology` replaces the gamma map with `unit_gammas()` for the resolved tree, which sets `gamma = 1.0` on every edge ([packages/treetime/src/timetree/inference/time_inference.rs](../../packages/treetime/src/timetree/inference/result.rs))
- The next round recomputes the multipliers at its start, so the reset lasts for the rest of the round that resolved polytomies, and for the final outputs when that round is the last one

## v0 behavior

- Same order within a round: relaxed clock, then polytomies ([packages/legacy/treetime/treetime/treetime.py#L311-L343](../../packages/legacy/treetime/treetime/treetime.py#L311-L343))
- After resolution, `make_time_tree` rebuilds the branch-length interpolators but copies each existing node's gamma first ([packages/legacy/treetime/treetime/clock_tree.py#L349-L370](../../packages/legacy/treetime/treetime/clock_tree.py#L349-L370)). Only new nodes start at 1.0 ([packages/legacy/treetime/treetime/branch_len_interpolator.py#L29](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L29))

## Impact

- With `--relax` and `--resolve-polytomies`, the relaxed clock has no effect in rounds that resolve polytomies
- When the last round resolves polytomies, the reported gamma of every edge is 1.0

## Open question

v0 keeps gamma on existing edges and sets 1.0 on new edges. The alternative is to refit the relaxed clock on the resolved tree before the rebuild. Either changes outputs and needs approval.
