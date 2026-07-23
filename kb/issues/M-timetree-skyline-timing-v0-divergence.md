# Timetree skyline coalescent timing may diverge from v0

v1 re-fits the coalescent prior (constant or skyline) once per optimization round on that round's node times, so the output node times are conditioned on a prior fit on refined times. After the loop and any final marginal pass, it re-fits the prior once more on the final node times and reports that value, so the reported prior is the maximum-likelihood fit of the tree actually written out. v0 activates the skyline on the last iteration within the loop (constant Tc on earlier iterations), so its final node times are conditioned on the skyline prior during the loop itself.

The v0 code (`packages/legacy/treetime/treetime/treetime.py` line 312) uses `if Tc == 'skyline' and niter < max_iter - 1: tmpTc = 'const'`, so the last iteration (`niter == max_iter - 1`) runs with `Tc = 'skyline'` inside the loop. v1 instead re-fits every round rather than switching from constant to skyline only on the last iteration.

Both approaches produce self-consistent results. The remaining difference is the constant-Tc-until-last-iteration schedule: v1 uses the skyline (or constant optimum) from round 1, v0 uses a constant Tc until the final iteration. For datasets where the skyline shape deviates from constant Tc, final node time estimates may differ.

## Investigation needed

- Run v0 and v1 on a dataset with population size variation (e.g. `flu/h3n2/200` with `--coalescent-skyline`)
- Compare node times to determine whether the difference is measurable
- If measurable, decide whether v1's approach is an intentional improvement or should match v0

## Related

- `kb/issues/N-coalescent-skyline-simplex-initialization-undecided.md`
- `kb/issues/N-coalescent-skyline-extrapolation-policy-undecided.md`
- `kb/issues/N-coalescent-skyline-quadrature-contract-undecided.md`
