# Hoist indel rate estimation before timetree refinement loop

## Summary

The timetree inference runner recomputes `estimate_indel_rate()` on every call to `run_timetree`, causing the indel rate to oscillate across refinement iterations as branch lengths change.

## Details

`packages/treetime/src/timetree/inference/runner.rs:95:`

Each refinement iteration calls `run_timetree` which calls `estimate_indel_rate(graph, partitions)`. Since branch lengths change between iterations (that is the purpose of refinement), the indel rate estimate changes too, creating a feedback loop: branch length changes shift the indel rate, which shifts the Poisson indel likelihood contribution, which shifts branch length targets.

The optimize command hoisted indel rate estimation before its loop (PR #619), but the timetree runner still recomputes per iteration.

## Impact

- Indel rate oscillates across timetree refinement iterations
- Convergence behavior is less predictable due to the feedback loop
- Branch length estimates in indel-rich regions are less stable

## Fix

Hoist `estimate_indel_rate()` before the timetree refinement loop, matching the optimize command pattern. Compute once and hold fixed across iterations.

## Related issues

- Source: [N-timetree-indel-rate-recomputed-per-iteration.md](../issues/N-timetree-indel-rate-recomputed-per-iteration.md) -- delete after full resolution
