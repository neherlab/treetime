# Node dates move by up to 0.15 years when a time grid gains or loses one point

A node date is the grid point with the highest posterior density: `fn DistributionFunction::likely_time()` takes the argmax of the stored ordinates [packages/treetime-distribution/src/distribution_core/function.rs#L297-L314](../../packages/treetime-distribution/src/distribution_core/function.rs#L297-L314) (the root takes the posterior median instead, [kb/decisions/timetree-root-time-at-posterior-median.md](../decisions/timetree-root-time-at-posterior-median.md)). For a weakly determined node the posterior is flat near its peak, and a small change of the grid moves the argmax far.

## Evidence

Commit `559c8697f` changed only the point count of the mass window in `fn resample_to_mass_window()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L332-L368](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L332-L368): a rounded-up quotient that is an integer in exact arithmetic gave 300 or 301 points, and the fix counts it exactly. Smoke runs before and after the fix (2026-10-19):

- **Most cases**: the Auspice `num_date` of every node moved by at most 0.008 years in 11 of the 13 timetree cases with changed dates
- **`zika/20` with `--coalescent-opt`**: one node moved by 0.152 years; its branch length changed by a relative 0.935, and the fitted $T_c$ moved from 5.4815 to 5.4824
- **`ebola/362` with the skyline coalescent**: dates moved by up to 0.072 years

A repeat of the `zika/20` reproduction with `--coalescent-opt` on the code of 2026-10-20 (mass-model posterior statistics, forward-posterior trim) moved node dates by at most 0.0023 years, against 90% widths of 0.25 to 0.84 years, so there the moves lie far inside the stated uncertainty.

> [!IMPORTANT]
> **Investigation required.** Whether the 0.152-year move of the original measurement stayed within that node's 90% interval or jumped between separate modes of its posterior. The repeat does not reproduce a move of that size, so the question needs the posterior of a node whose date moves by a large part of its width under a one-point grid change.

## Fix direction

> [!IMPORTANT]
> **Decision required.** A date that does not depend on single grid points. v0 also dates a node at the grid point of its highest posterior density ([packages/legacy/treetime/treetime/clock_tree.py#L922](../../packages/legacy/treetime/treetime/clock_tree.py#L922), [packages/legacy/treetime/treetime/distribution.py#L210-L212](../../packages/legacy/treetime/treetime/distribution.py#L210-L212)), so every option except (c) diverges from v0. Options: (a) the peak of a local fit around the grid argmax, which keeps the meaning of the date as the posterior mode and moves it by less than one grid cell; (b) the posterior median for every node, as the root already uses, which changes the meaning of the date; (c) no change, with the sensitivity stated in the documentation of the date outputs.

## Related

- [H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md](H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md)
