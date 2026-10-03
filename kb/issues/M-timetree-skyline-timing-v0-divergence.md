# Timetree skyline coalescent timing diverges from v0

v0 and v1 apply the skyline coalescent prior at different points of the refinement loop, and only v0 disables the early convergence exit in skyline mode. No decision records either difference.

## v0 behavior

- Inside the loop, every iteration except the last uses a constant $T_c$ prior; only the last iteration (`niter == max_iter - 1`) fits and applies the skyline: `if Tc == 'skyline' and niter < max_iter - 1: tmpTc = 'const'` ([packages/legacy/treetime/treetime/treetime.py#L312-L315](../../packages/legacy/treetime/treetime/treetime.py#L312-L315))
- The loop never stops early in skyline mode: the convergence check requires `Tc != 'skyline'` ([packages/legacy/treetime/treetime/treetime.py#L374](../../packages/legacy/treetime/treetime/treetime.py#L374)), so a skyline run always performs `max_iter` iterations

## v1 behavior

- The initial round fits the skyline from the first node times (`fn setup_coalescent` in [packages/treetime/src/timetree/round.rs](../../packages/treetime/src/timetree/round.rs)), and every loop round re-fits it from the current node times before the round runs ([packages/treetime/src/timetree/refinement_loop.rs#L54-L63](../../packages/treetime/src/timetree/refinement_loop.rs#L54-L63)), so every round uses a skyline prior
- No skyline step runs after the loop
- The early convergence exit applies in skyline mode as in every other mode (`fn TimetreeOptimizer::next_iter` in [packages/treetime/src/timetree/convergence/optimizer.rs](../../packages/treetime/src/timetree/convergence/optimizer.rs)), so a skyline run can stop before `max_iter` with a skyline fitted to times that v0 would have refined further

## Impact

For datasets where the skyline shape deviates from a constant $T_c$, the intermediate node times differ from v0 because each v1 round already conditions on the skyline. A skyline run that converges early performs fewer iterations than v0, so final node times and the fitted skyline may differ.

## Investigation needed

- Run v0 and v1 on a dataset with population size variation (e.g. `flu/h3n2/200` with `--coalescent-skyline`)
- Compare node times, the fitted skyline, and the iteration count to determine whether the differences are measurable
- If measurable, decide whether v1's schedule is an intentional change or should match v0

## Related

- `kb/decisions/coalescent-skyline-convex-log-tc.md`
- `kb/issues/N-coalescent-skyline-extrapolation-policy-undecided.md`
- `kb/issues/N-coalescent-skyline-quadrature-contract-undecided.md`
