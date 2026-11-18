# Timetree indel rate recomputed per refinement iteration

## Summary

Every call to `run_timetree` re-estimates the per-partition indel rates from the current branch lengths, so the indel rates change across refinement iterations as branch lengths change.

## Details

With marginal branch lengths, `fn run_timetree()` reaches `fn compute_branch_distributions_marginal_mode()`, which calls `IndelPolicy::edge_indels(graph, partitions, branch_lengths)` ([packages/treetime/src/timetree/inference/runner.rs#L329](../../packages/treetime/src/timetree/inference/runner.rs#L329)). With `IndelPolicy::Include`, `edge_indels()` calls the private `fn estimate_indel_rate()` once per partition: total indels divided by total branch length ([packages/treetime/src/optimize/indel.rs#L24-L44](../../packages/treetime/src/optimize/indel.rs#L24-L44), [packages/treetime/src/optimize/indel.rs#L148-L168](../../packages/treetime/src/optimize/indel.rs#L148-L168)).

The callers of `run_timetree` are:

- `fn initialize()` ([packages/treetime/src/timetree/init.rs#L75](../../packages/treetime/src/timetree/init.rs#L75))
- every refinement round, through `fn refresh_times()` and `fn infer_times()` ([packages/treetime/src/timetree/round.rs#L255](../../packages/treetime/src/timetree/round.rs#L255))
- `fn infer_final_times()` ([packages/treetime/src/timetree/round.rs#L82](../../packages/treetime/src/timetree/round.rs#L82))
- the confidence reruns of `fn compute_rate_susceptibility()` ([packages/treetime/src/timetree/confidence.rs#L51](../../packages/treetime/src/timetree/confidence.rs#L51))

Since branch lengths change between iterations (that is the purpose of refinement), the indel rate estimate changes too, creating a feedback loop: branch length changes shift the indel rate, which shifts the Poisson indel likelihood contribution, which shifts branch length targets.

The optimize command estimates the rates once before its loop ([packages/treetime/src/optimize/run_loop.rs#L47](../../packages/treetime/src/optimize/run_loop.rs#L47)); `EdgeIndels::recount()` then refreshes the counts and keeps the rates ([packages/treetime/src/optimize/indel.rs#L57-L65](../../packages/treetime/src/optimize/indel.rs#L57-L65)). The timetree runner has no such split.

## Impact

- Indel rate oscillates across timetree refinement iterations
- Convergence behavior is less predictable due to the feedback loop
- Branch length estimates in indel-rich regions are less stable

## Fix

Estimate the indel rates once before the timetree refinement loop, matching the optimize command, and pass the `EdgeIndels` into the time inference. Hold the rates fixed across iterations and recount the indels with `EdgeIndels::recount()` when the reconstruction changes.
