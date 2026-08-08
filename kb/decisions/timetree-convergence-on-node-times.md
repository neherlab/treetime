# Timetree convergence is judged on node-time movement

The timetree refinement loop converges when no node time moved more than
`NODE_TIME_TOLERANCE_YEARS = 1e-2` during a round and no polytomy was resolved. The previous
criterion, `n_diff == 0 && n_resolved == 0`, counted changes in the reconstructed ancestral
sequences, which is not what the loop moves.

**Type**: New v1 criterion. v0's timetree loop uses `ndiff == 0 && n_resolved == 0`.

**v1 location**:
[`ConvergenceMetrics::has_converged`](../../packages/treetime/src/timetree/convergence/metrics.rs),
measured in [`Refinement::run`](../../packages/treetime/src/timetree/refinement.rs) via
[`convergence/node_times.rs`](../../packages/treetime/src/timetree/convergence/node_times.rs).

## Why `n_diff` cannot see the loop

`n_diff` counts positions where the maximum-likelihood ancestral state differs between rounds. That
state rarely flips even when every node time shifts. Measured on `data/ebola/20` after clock-
constrained propagation was in place: node dates moved by up to 0.023 y in round 1 and 0.065 y in
round 2, and `n_diff` was 0 in every round. Under the old criterion the loop stopped after round 1
on all three configurations tested while the tree was still visibly moving.

The criterion is therefore node-time movement, with `n_diff` retained only as the fallback for a
tree where no node is dated on both sides of the round and no movement is measurable. `n_resolved
== 0` remains required: a round that changed the topology is not converged whatever the times did.

## What is measured

Times are snapshotted at the top of `Refinement::run` and compared after `rebuild_inference`, so
the quantity is the movement *within* a round rather than a difference between rounds. This is
well-defined from round 1 and needs no seeding. Nodes present in only one snapshot are skipped;
polytomy resolution introduces nodes with no earlier position, and such a round is excluded by
`n_resolved` anyway.

Both `max_time_change` and `rms_time_change` are recorded and written to the tracelog CSV. Their
ratio is diagnostic: `max` far above `rms` means one node is oscillating while the rest of the tree
is still, i.e. local bistability rather than a loop that failed to settle. A ratio of
$\sqrt{n_{\text{nodes}}}$ is exactly one node moving.

## The tolerance

`1e-2` years is about 3.7 days. It is set at the resolution the data carry, not the resolution the
solver can report:

- sampling dates are rarely recorded finer than a day;
- the biological replication cycle is itself a few days, so movement below that is not
  interpretable;
- node times are quantized by the 300-point branch-length grid. Movements come out as integer
  multiples of a per-dataset step — about 0.0026 y on `data/ebola/20`, about 0.016 y on
  `data/mpox/clade-ii/1000`. A tolerance below the step makes runs exhaust `--max-iter` chasing
  jitter that no longer changes the answer.

A tolerance derived from the grid step rather than fixed is the obvious refinement; see
[N-timetree-convergence-tolerance-vs-branch-grid.md](../issues/N-timetree-convergence-tolerance-vs-branch-grid.md).

## Results

`data/ebola/20`, no coalescent, converging in 3 rounds where the old criterion stopped after 1:

```
Iteration 1: max_dt=0.0229, rms_dt=0.0046, n_diff=0, n_resolved=0
Iteration 2: max_dt=0.0652, rms_dt=0.0329, n_diff=0, n_resolved=0
Iteration 3: max_dt=0.0064, rms_dt=0.0018, n_diff=0, n_resolved=0  [converged]
```

## Related

- [timetree-clock-constrained-profile-propagation.md](timetree-clock-constrained-profile-propagation.md)
  — without it there is no movement to measure
- [M-timetree-convergence-metric-deficiencies.md](../issues/M-timetree-convergence-metric-deficiencies.md)
  — the remaining defects in the metric set
