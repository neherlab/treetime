# Timetree skyline coalescent timing diverges from v0

v0 and v1 apply the skyline coalescent prior at different points of the refinement loop, and only v0 disables the early convergence exit in skyline mode. No decision records either difference.

## v0 behavior

- Inside the loop, every iteration except the last uses a constant $T_c$ prior; only the last iteration (`niter == max_iter - 1`) fits and applies the skyline: `if Tc == 'skyline' and niter < max_iter - 1: tmpTc = 'const'` ([packages/legacy/treetime/treetime/treetime.py#L312-L315](../../packages/legacy/treetime/treetime/treetime.py#L312-L315))
- The loop never stops early in skyline mode: the convergence check requires `Tc != 'skyline'` ([packages/legacy/treetime/treetime/treetime.py#L374](../../packages/legacy/treetime/treetime/treetime.py#L374)), so a skyline run always performs `max_iter` iterations

## v1 behavior

- The init phase fits the skyline to the expected lineage counts of the Gaussian node marginals (`fn initialize()` [packages/treetime/src/timetree/init.rs#L84](../../packages/treetime/src/timetree/init.rs#L84)), and every loop round re-fits it before the round runs ([packages/treetime/src/timetree/refinement_loop.rs#L47-L49](../../packages/treetime/src/timetree/refinement_loop.rs#L47-L49)), from the expected lineage counts of the posteriors of the previous round ([packages/treetime/src/timetree/refinement_loop.rs#L71](../../packages/treetime/src/timetree/refinement_loop.rs#L71)), so every round uses a skyline prior
- No skyline step runs after the loop
- The early convergence exit applies in skyline mode as in every other mode (`fn TimetreeOptimizer::next_iter()` [packages/treetime/src/timetree/convergence/optimizer.rs#L35](../../packages/treetime/src/timetree/convergence/optimizer.rs#L35)), so a skyline run can stop before `max_iter` with a skyline fitted to times that v0 would have refined further

- Smoothing strength: the skyline objective adds a penalty on $z_i = \ln T_{c,i}$ with stiffness $s$ (default 2 in both versions), $\frac{s}{2} \sum_i (z_{i+1} - z_i)^2$ in v1 (`fn skyline_cost()` in `packages/treetime/src/coalescent/skyline.rs`) and $s \sum_i (z_{i+1} - z_i)^2$ in v0 ([packages/legacy/treetime/treetime/merger_models.py#L312](../../packages/legacy/treetime/treetime/merger_models.py#L312)), so v1 smooths half as strongly at the same default. [kb/decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md) does not record the factor. On simulated data with known $T_c$ the v1 skyline is nearly unbiased (`keep-skyline-fixed`, $L = 10^4$, 30 runs on 2026-10-24: mean log bias +0.023, band coverage 0.97)
- Shape of $T_c$: v1 is piecewise constant on equal-width segments; v0 interpolates between evenly spaced points ([packages/legacy/treetime/treetime/merger_models.py#L305](../../packages/legacy/treetime/treetime/merger_models.py#L305))

## Impact

For datasets where the skyline shape deviates from a constant $T_c$, the intermediate node times differ from v0 because each v1 round already conditions on the skyline. A skyline run that converges early performs fewer iterations than v0, so final node times and the fitted skyline may differ.

## Investigation needed

- Run v0 and v1 on a dataset with population size variation (e.g. `flu/h3n2/200` with `--coalescent skyline`)
- Compare node times, the fitted skyline, and the iteration count to determine whether the differences are measurable
- If measurable, decide whether v1's schedule is an intentional change or should match v0

## Related

- `kb/decisions/coalescent-skyline-convex-log-tc.md`
- `kb/issues/N-coalescent-skyline-extrapolation-policy-undecided.md`
- `kb/issues/N-coalescent-skyline-quadrature-contract-undecided.md`
