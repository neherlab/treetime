# Timetree skyline coalescent timing may diverge from v0

v1 runs skyline coalescent optimization after the iteration loop exits, then re-runs timetree with the fitted skyline prior. v0 activates skyline on the last iteration within the loop, so the timetree inference runs once with the skyline prior during the loop itself.

The v1 code at `packages/treetime/src/timetree/pipeline.rs` (skyline block after the optimizer loop) has a comment claiming "matching v0 behavior where skyline optimization happens only on the final iteration after times have converged." The v0 code (`packages/legacy/treetime/treetime/treetime.py` line 312) uses `if Tc == 'skyline' and niter < max_iter - 1: tmpTc = 'const'`, which means the last iteration (when `niter == max_iter - 1`) runs with `Tc = 'skyline'` inside the loop, not after it. The comment may be inaccurate.

Both approaches produce self-consistent results. The difference is whether the final node times are conditioned on the skyline prior during the optimization iteration or via a separate post-loop pass. For datasets where the skyline shape deviates from constant Tc, final node time estimates may differ.

## Investigation needed

- Run v0 and v1 on a dataset with population size variation (e.g. `flu/h3n2/200` with `--coalescent-skyline`)
- Compare node times to determine whether the difference is measurable
- If measurable, decide whether v1's approach is an intentional improvement or should match v0

## Related

- `kb/issues/N-coalescent-skyline-simplex-initialization-undecided.md`
- `kb/issues/N-coalescent-skyline-extrapolation-policy-undecided.md`
- `kb/issues/N-coalescent-skyline-quadrature-contract-undecided.md`
