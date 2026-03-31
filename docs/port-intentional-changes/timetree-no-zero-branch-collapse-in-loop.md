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

Topology cleanup in the substitution-space M-step is consistent with EM algorithm theory <a id="cite-1"></a>[Dempster, Laird, and Rubin 1977](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x) [[1](#ref-1)]: each M-step maximizes the expected complete-data log-likelihood under its own parameterization. The substitution-space M-step optimizes branch lengths freely and can detect zero-length branches. The time-space M-step optimizes node times under clock constraints, where zero time difference is a valid state.

## Cross-references

- [Command relationships](../reports/command-relationships/_index.md): Documents the two sibling loops (optimize and timetree) sharing an E-step.

## References

1. <a id="ref-1"></a> Dempster, Arthur P., Nan M. Laird, and Donald B. Rubin. 1977. "Maximum Likelihood from Incomplete Data Via the EM Algorithm." _Journal of the Royal Statistical Society: Series B (Methodological)_ 39(1):1-22. https://doi.org/10.1111/j.2517-6161.1977.tb01600.x [↩](#cite-1)
