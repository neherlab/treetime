# Relaxed clock and polytomy mutation counts read the unconstrained branch length

Timetree edges carry two branch lengths: `branch_length`, the free ML or input estimate, and
`clock_branch_length`, the length the inferred node times imply. See
[timetree-clock-constrained-profile-propagation.md](../decisions/timetree-clock-constrained-profile-propagation.md).
Two consumers that want the clock-constrained value still read the free one.

## apply_relaxed_clock has almost no signal to fit

`fn apply_relaxed_clock()`
[packages/treetime/src/timetree/optimization/relaxed_clock.rs#L27-L36](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs#L27-L36)
compares `act_len = time_length * clock_rate` against `opt_len = branch_length`. Both are
ML-derived: `time_length` is set from `distribution.likely_time()`, the unconstrained per-edge
peak. So the ratio the per-branch rate multiplier $\gamma$ is fitted from is close to 1 by
construction and carries little information about clock rate variation.

v0 compares `clock_length` against `mutation_length`, i.e. the constrained length against the free
one, which is the contrast $\gamma$ is meant to capture. The v1 equivalent is
`clock_branch_length` against `branch_length`. `struct RoundState` already holds the constrained lengths in
`clock_branch_lengths`
([packages/treetime/src/timetree/round.rs#L346](../../packages/treetime/src/timetree/round.rs#L346)).

## Polytomy mutation counts fall back to the free length

`fn edge_mutation_count()`
[packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L192-L216](../../packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L192-L216)
estimates a branch's substitution count as `round(branch_length * total_length)` when the
reconstructed substitution list cannot be read. The caller passes the free `branch_length`
([packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L175-L180](../../packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L175-L180)).
The estimate should use the clock-constrained length for consistency with the sweep it feeds,
which places those substitutions in time against a coalescent rate defined on the same time axis.

The fallback runs on every edge in input branch-length mode (`BranchModel::Input` has no
reconstruction), and in marginal mode on an edge for which a sparse partition stores no
substitution estimates.

## Not affected

The transition counts of GTR inference (`fn count_transitions_dense()` in
[packages/treetime/src/partition/marginal/shared/data.rs#L15](../../packages/treetime/src/partition/marginal/shared/data.rs#L15)
and `fn count_transitions_sparse()` in
[packages/treetime/src/partition/marginal/sparse/count.rs#L16](../../packages/treetime/src/partition/marginal/sparse/count.rs#L16))
also read the ML length. That is unobservable: GTR inference runs only before the refinement loop,
so the two lengths are equal at every point where the counts are used. It becomes a real choice only
if GTR re-estimation moves into the loop.

## Impact

- Relaxed clock (`--relax`) fits $\gamma$ against a contrast that is near-degenerate, so the
  feature does substantially less than intended. Not quantified.
- Polytomy resolution mis-weights branches in the fallback path, which input branch-length mode
  takes on every edge. Input mode fails earlier for an unrelated reason, see
  [H-timetree-input-branch-lengths-abort-on-point-division.md](H-timetree-input-branch-lengths-abort-on-point-division.md).

## Validation

The relaxed clock change alters $\gamma$ on every edge and therefore node times; it needs a
before/after comparison on `data/ebola/20` and `data/mpox/clade-ii/1000` rather than a unit test
alone.
