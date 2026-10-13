# The optimize loop returns the last tree when the likelihood declines from the start

## Problem

`fn run_optimize_loop` [packages/treetime/src/optimize/run_loop.rs#L32](../../packages/treetime/src/optimize/run_loop.rs#L32) records the best branch lengths [packages/treetime/src/optimize/run_loop.rs#L90-L93](../../packages/treetime/src/optimize/run_loop.rs#L90-L93), but reverts to them only when a drop follows a new best [packages/treetime/src/optimize/run_loop.rs#L105](../../packages/treetime/src/optimize/run_loop.rs#L105):

```rust
if i >= 2 && iteration_lh.total_lh < lh_prev && lh_prev >= best_lh {
```

When the first iteration is the best and the likelihood declines after it, the condition never holds. The loop runs to `max_iter` and returns the last branch lengths, which have the lowest likelihood.

## Evidence

`optimize` on simulated data of the validation suite ([validation/README.md](../../validation/README.md), schedule `exp`, replicate `rep07`, alignment length $L = 1000$, substitution rate $\mu = 0.01$ per site and year, started from branch lengths with a higher likelihood than its result) lowered the log likelihood at every iteration from -11542.8 to -11565.5 and returned the final tree. This run used the branch objective before its transition matrix was oriented as child given parent. Simultaneous per-branch updates need not increase the likelihood, so the same pattern can occur with the corrected objective.

## Fix direction

Return the best branch lengths seen whenever the loop ends.

## Validation

A regression test with a likelihood that declines from the first iteration: the loop returns the branch lengths of the first iteration.
