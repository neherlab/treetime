# Clock branch-length damping is a flat factor, unlike the optimize loop's schedule

`CLOCK_BRANCH_LENGTH_DAMPING = 0.5` in
[`packages/treetime/src/timetree/inference/runner.rs`](../../packages/treetime/src/timetree/inference/runner.rs)
is applied unchanged in every refinement round. The branch-length optimization loop instead uses an
exponential schedule in
[`apply_damping`](../../packages/treetime/src/optimize/iteration.rs): `f = max(d^(i+1), 0.01)` with
`d = 0.75`, so early rounds take conservative steps and late rounds approach the full update.

A flat factor never approaches the undamped step, so the final rounds converge geometrically at
rate `1 - f` rather than reaching the fixed point directly. The fixed point is the same either way
($b = (1-f)b + fb$ for any $f$), so this is a rate question, not a correctness one.

## Why it is flat today

`commit_clock_branch_lengths` is called from four sites, only one of which is inside the loop. The
round index is not in scope at the commit and would have to be threaded through `Refinement` to
build a schedule. That plumbing was not worth adding before knowing whether the flat factor costs
measurable rounds.

## Evidence

`data/ebola/20` converges in 3 rounds with no coalescent and 5 with `--coalescent-opt`, well inside
the default `--max-iter`. No dataset has yet been shown to need more rounds because of the flat
factor.

## Impact

None demonstrated. Possible extra rounds on datasets that converge slowly.

## Options

- Thread the round index into the commit and use the same `d^(i+1)` schedule as the optimize loop,
  for consistency between the two loops.
- Leave flat and revisit if a dataset is found whose round count is damping-limited, which would
  show as a `max_time_change` decaying by a constant ratio near 0.5 per round.
