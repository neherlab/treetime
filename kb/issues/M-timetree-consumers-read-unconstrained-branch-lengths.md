# Relaxed clock and polytomy mutation counts read the unconstrained branch length

Timetree edges now carry two branch lengths: `branch_length`, the free ML or input estimate, and
`clock_branch_length`, the length the inferred node times imply. See
[timetree-clock-constrained-profile-propagation.md](../decisions/timetree-clock-constrained-profile-propagation.md).
Two consumers that want the clock-constrained value still read the free one.

## apply_relaxed_clock has almost no signal to fit

[`packages/treetime/src/timetree/optimization/relaxed_clock.rs`](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs)
compares `act_len = time_length * clock_rate` against `opt_len = branch_length`. Both are
ML-derived — `time_length` is set from `distribution.likely_time()`, the unconstrained per-edge
peak — so the ratio the per-branch rate multiplier $\gamma$ is fitted from is close to 1 by
construction and carries little information about clock rate variation.

v0 compares `clock_length` against `mutation_length`, i.e. the constrained length against the free
one, which is the contrast $\gamma$ is meant to capture. The v1 equivalent is
`clock_branch_length` against `branch_length`.

## Polytomy mutation counts fall back to the free length

`edge_mutation_count` in
[`packages/treetime/src/timetree/optimization/polytomy/mod.rs`](../../packages/treetime/src/timetree/optimization/polytomy/mod.rs)
estimates a branch's substitution count as `round(branch_length * total_length)` when the
reconstructed substitution list cannot be read. The estimate should use the clock-constrained
length for consistency with the sweep it feeds, which places those substitutions in time against a
coalescent rate defined on the same time axis.

## Not affected

`count_transitions` in
[`marginal_core.rs`](../../packages/treetime/src/partition/marginal_core.rs) and
[`marginal_sparse.rs`](../../packages/treetime/src/partition/marginal_sparse.rs) also reads the ML
length. That is currently unobservable: GTR inference only runs before the refinement loop, so the
two lengths are equal at every point where the counts are used. It becomes a real choice only if
GTR re-estimation moves into the loop.

## Impact

- Relaxed clock (`--relax`) fits $\gamma$ against a contrast that is near-degenerate, so the
  feature does substantially less than intended. Not quantified.
- Polytomy resolution mis-weights branches in the fallback path only, which is reached when the
  marginal forward pass has not yet run.

## Validation

The relaxed clock change alters $\gamma$ on every edge and therefore node times; it needs a
before/after comparison on `data/ebola/20` and `data/mpox/clade-ii/1000` rather than a unit test
alone.
