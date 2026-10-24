# Node dates move by up to 0.15 years when a time grid gains or loses one point

A node date is the grid point with the highest posterior density: `fn DistributionFunction::likely_time()` takes the argmax of the stored ordinates [packages/treetime-distribution/src/distribution_core/function.rs#L297-L314](../../packages/treetime-distribution/src/distribution_core/function.rs#L297-L314) (the root takes the posterior median instead, [kb/decisions/timetree-root-time-at-posterior-median.md](../decisions/timetree-root-time-at-posterior-median.md)). For a weakly determined node the posterior is flat near its peak, and a small change of the grid moves the argmax far.

## Evidence

Commit `559c8697f` changed only the point count of the mass window in `fn resample_to_mass_window()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L332-L368](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L332-L368): a rounded-up quotient that is an integer in exact arithmetic gave 300 or 301 points, and the fix counts it exactly. Smoke runs before and after the fix (2026-10-19):

- **Most cases**: the Auspice `num_date` of every node moved by at most 0.008 years in 11 of the 13 timetree cases with changed dates
- **`zika/20` with `--coalescent opt`**: one node moved by 0.152 years; its branch length changed by a relative 0.935, and the fitted $T_c$ moved from 5.4815 to 5.4824
- **`ebola/362` with the skyline coalescent**: dates moved by up to 0.072 years

A repeat of the `zika/20` reproduction with `--coalescent opt` on the code of 2026-10-20 (mass-model posterior statistics, forward-posterior trim) moved node dates by at most 0.0023 years, against 90% widths of 0.25 to 0.84 years, so there the moves lie far inside the stated uncertainty.

The same sensitivity appears without any change of a point count. A change of the mass model that alters cell masses only by a relative $u / \Delta x$ (each cell takes the width between its rounded end coordinates instead of the stored spacing, with $u$ the unit in the last place of the coordinate) moved dates on `data/sc2/5100` by up to 0.042 of a node's 90% width, and 1491 of 10 199 nodes by more than 0.1% of their width. The first difference is a one-ulp move of the parent of the leaf `USA/MN-MDH-48906/2025` in the first round; the later rounds of the refinement loop spread it. The same change moved dates of `data/ebola/362`, `data/flu/h3n2/200` and `data/zika/86` by at most $2.3 \cdot 10^{-10}$ of a width. That branch goes from $-1$ ulp to 0, which changes the count of the warning about branches that run backward in time (`fn blended_clock_branch_lengths()`, [packages/treetime/src/timetree/inference/runner.rs#L189-L243](../../packages/treetime/src/timetree/inference/runner.rs#L189-L243)) but not the committed branch length, which is 0 in both cases.

> [!IMPORTANT]
> **Investigation required.** Which discrete step turns the one-ulp difference on `data/sc2/5100` into moves of 0.042 of a width: the argmax of `fn DistributionFunction::likely_time()`, a rounded point count, a trim cut that snaps to a grid point, or a refinement round that converges to another state. Compare the per-round node dates of the two builds and name the first step whose output differs by more than rounding.

> [!IMPORTANT]
> **Investigation required.** Whether the 0.152-year move of the original measurement stayed within that node's 90% interval or jumped between separate modes of its posterior. The repeat does not reproduce a move of that size, so the question needs the posterior of a node whose date moves by a large part of its width under a one-point grid change.

## Fix direction

> [!IMPORTANT]
> **Decision required.** A date that does not depend on single grid points. v0 also dates a node at the grid point of its highest posterior density ([packages/legacy/treetime/treetime/clock_tree.py#L922](../../packages/legacy/treetime/treetime/clock_tree.py#L922), [packages/legacy/treetime/treetime/distribution.py#L210-L212](../../packages/legacy/treetime/treetime/distribution.py#L210-L212)), so every option except (c) diverges from v0. Options: (a) the peak of a local fit around the grid argmax, which keeps the meaning of the date as the posterior mode and moves it by less than one grid cell; (b) the posterior median for every node, as the root already uses, which changes the meaning of the date; (c) no change, with the sensitivity stated in the documentation of the date outputs.

## Related

- [discretization error of the default grid](../decisions/distribution-tails-and-arithmetic.md#discretization-error-of-the-default-grid)
