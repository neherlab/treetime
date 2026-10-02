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

Four places in `timetree/round.rs` compute clock branch lengths with `fn blended_clock_branch_lengths`. Three go through `fn RoundState::blend_clock_branch_lengths` [packages/treetime/src/timetree/round.rs#L191-L205](../../packages/treetime/src/timetree/round.rs#L191-L205):

- `fn run_initial_round` [packages/treetime/src/timetree/round.rs#L106](../../packages/treetime/src/timetree/round.rs#L106): undamped, before the refinement loop
- `fn refresh_times` [packages/treetime/src/timetree/round.rs#L364](../../packages/treetime/src/timetree/round.rs#L364): damped with `CLOCK_BRANCH_LENGTH_DAMPING`, inside the refinement loop
- `fn final_marginal_round` [packages/treetime/src/timetree/round.rs#L160](../../packages/treetime/src/timetree/round.rs#L160): undamped, after the refinement loop

`fn refine_topology` [packages/treetime/src/timetree/round.rs#L325-L333](../../packages/treetime/src/timetree/round.rs#L325-L333) calls it directly: undamped, inside the refinement loop, after a polytomy resolution changes the topology.

Only the `refresh_times` site applies the flat factor. `fn refinement_round` [packages/treetime/src/timetree/round.rs#L115-L122](../../packages/treetime/src/timetree/round.rs#L115-L122) does not receive the round index, so a schedule needs the index threaded through `refinement_round` and `refresh_times`. This plumbing waits until a dataset shows that the flat factor costs measurable rounds.

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
