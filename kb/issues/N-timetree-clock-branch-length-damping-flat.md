# Clock branch-length damping is a flat factor, unlike the optimize loop's schedule

`CLOCK_BRANCH_LENGTH_DAMPING = 0.5` [packages/treetime/src/timetree/inference/runner.rs#L183](../../packages/treetime/src/timetree/inference/runner.rs#L183) is applied unchanged in every refinement round: `fn blended_clock_branch_lengths()` sets each clock branch length to `(1 - damping) * previous + damping * fresh` [packages/treetime/src/timetree/inference/runner.rs#L209-L212](../../packages/treetime/src/timetree/inference/runner.rs#L209-L212). The branch-length optimization loop instead uses a schedule in `fn apply_damping()` [packages/treetime/src/optimize/iteration.rs#L7-L23](../../packages/treetime/src/optimize/iteration.rs#L7-L23): the old length gets the weight `max(d^(i+1), 0.01)` in round `i`, so early rounds take conservative steps and late rounds approach the full update. `d` is the `--damping` value of `treetime optimize`, default `0.0` (no damping).

The two parameters have opposite meanings: in `blended_clock_branch_lengths()` the damping value weights the fresh length (`CLOCK_BRANCH_LENGTH_UNDAMPED = 1.0` keeps only the fresh length), in `apply_damping()` the damping factor weights the old length.

A flat factor never approaches the undamped step, so the final rounds converge geometrically at rate 0.5 rather than reaching the fixed point directly. The fixed point is the same either way ($b = (1-f)b + fb$ for any $f$), so this is a rate question, not a correctness one.

## Why it is flat today

Four places in `timetree/` compute clock branch lengths with `fn blended_clock_branch_lengths()`. Three go through `fn RoundState::blend_clock_branch_lengths()` [packages/treetime/src/timetree/round.rs#L372-L386](../../packages/treetime/src/timetree/round.rs#L372-L386):

- `fn initialize()` [packages/treetime/src/timetree/init.rs#L100](../../packages/treetime/src/timetree/init.rs#L100): undamped, before the refinement loop
- `fn refresh_times()` [packages/treetime/src/timetree/round.rs#L215-L235](../../packages/treetime/src/timetree/round.rs#L215-L235): damped with `CLOCK_BRANCH_LENGTH_DAMPING`, inside the refinement loop
- `fn final_marginal_round()` [packages/treetime/src/timetree/round.rs#L88-L99](../../packages/treetime/src/timetree/round.rs#L88-L99): undamped, after the refinement loop

`fn refine_topology()` [packages/treetime/src/timetree/round.rs#L126](../../packages/treetime/src/timetree/round.rs#L126) calls it directly [packages/treetime/src/timetree/round.rs#L175-L183](../../packages/treetime/src/timetree/round.rs#L175-L183): undamped, inside the refinement loop, after a polytomy resolution changes the topology.

Only the `refresh_times()` site applies the flat factor. `fn refinement_round()` [packages/treetime/src/timetree/round.rs#L39-L47](../../packages/treetime/src/timetree/round.rs#L39-L47) does not receive the round index, so a schedule needs the index passed through `refinement_round()` and `refresh_times()`.

## Evidence

`data/ebola/20` converges in 3 rounds with no coalescent and 5 with `--coalescent opt`, well inside the default `--max-iter`.

> [!IMPORTANT]
> **Investigation required.** No dataset has been shown to need more rounds because of the flat factor. Evidence to collect: the per-round `max_time_change` of the tracelog (`--output-tracelog`) on slowly converging datasets, for example `data/mpox/clade-ii/1000` and `data/sc2/4500`. A round count limited by the damping shows as a `max_time_change` that decays by a constant ratio near 0.5 per round.

## Impact

None demonstrated. Possible extra rounds on datasets that converge slowly.

## Options

- Pass the round index to the blend and use a decaying schedule like the optimize loop's, for consistency between the two loops
- Leave flat and revisit if the investigation above finds a dataset whose round count is damping-limited
