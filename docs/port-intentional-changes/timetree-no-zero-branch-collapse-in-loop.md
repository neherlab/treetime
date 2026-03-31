# Timetree EM loop does not collapse zero-length branches

## Deviation

The v1 timetree EM loop performs temporal polytomy resolution but does not collapse zero-length branches or merge shared mutations inside the loop. The optimize command's `prune_and_merge_in_loop()` is not called during timetree iterations.

## v0 behavior

v0 also does not prune inside the timetree EM loop. The `_run()` method at `packages/legacy/treetime/treetime/treetime.py:307:` iterates `make_time_tree()` + `infer_ancestral_sequences()` with polytomy resolution, but does not call `prune_short_branches()`. Pruning happens in the pre-loop `optimize_tree(max_iter=1)` calls (lines 243, 266), not during timetree iterations.

Inside the EM loop, when polytomies are resolved (line 329), v0 explicitly sets `prune_short=False` before calling `optimize_tree(max_iter=0)`.

## Rationale

The timetree M-step operates in calendar-time space with clock constraints ($b_e = \Delta t \cdot \mu \cdot \gamma$). Branch lengths are derived from node time differences, not optimized freely. Zero-length branches in the timetree context represent contemporaneous ancestors, which is a valid biological state (e.g. a common ancestor sampled at the same time as its descendants). Collapsing these would remove legitimate topology.

The optimize M-step operates in substitution space with free branch lengths ($b_e \geq 0$). A zero-length branch means no evolutionary distance, which is the condition for collapse.

The two loops are siblings sharing an E-step (ancestral reconstruction via `update_marginal()`) but with different M-step objectives. Topology cleanup belongs in the substitution-space objective, not the time-space objective.

Initial branch quality is handled by the optimize pre-step before the timetree loop begins (tracked in [M-timetree-missing-initial-branch-optimization](../port-known-issues/M-timetree-missing-initial-branch-optimization.md)).
