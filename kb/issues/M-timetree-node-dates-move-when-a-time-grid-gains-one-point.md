# Node dates move by up to 0.15 years when a time grid gains or loses one point

A node date is the grid point with the highest posterior density: `fn DistributionFunction::likely_time()` takes the argmax of the stored ordinates [packages/treetime-distribution/src/distribution_core/function.rs#L298-L315](../../packages/treetime-distribution/src/distribution_core/function.rs#L298-L315) (the root takes the posterior median instead, [kb/decisions/timetree-root-time-at-posterior-median.md](../decisions/timetree-root-time-at-posterior-median.md)). For a weakly determined node the posterior is flat near its peak, and a small change of the grid moves the argmax far.

## Evidence

Commit `559c8697f` changed only the point count of the mass window in `fn resample_to_mass_window()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L115-L121](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L115-L121): a rounded-up quotient that is an integer in exact arithmetic gave 300 or 301 points, and the fix counts it exactly. Smoke runs before and after the fix (2026-10-19):

- **Most cases**: the Auspice `num_date` of every node moved by at most 0.008 years in 11 of the 13 timetree cases with changed dates
- **`zika/20` with `--coalescent-opt`**: one node moved by 0.152 years; its branch length changed by a relative 0.935, and the fitted $T_c$ moved from 5.4815 to 5.4824
- **`ebola/362` with the skyline coalescent**: dates moved by up to 0.072 years

> [!IMPORTANT]
> **Investigation required.** Whether the moved nodes have flat posteriors whose 90% interval is much wider than the move, so the change is within the stated uncertainty, or whether the argmax jumps between separate modes. Compare the posterior of the moved `zika/20` node before and after the fix, and its 90% interval.

## Fix direction

> [!IMPORTANT]
> **Decision required.** A date that does not depend on single grid points. v0 also dates a node at the grid point of its highest posterior density ([packages/legacy/treetime/treetime/clock_tree.py#L922](../../packages/legacy/treetime/treetime/clock_tree.py#L922), [packages/legacy/treetime/treetime/distribution.py#L210-L212](../../packages/legacy/treetime/treetime/distribution.py#L210-L212)), so every option except (c) diverges from v0. Options: (a) the peak of a local fit around the grid argmax, which keeps the meaning of the date as the posterior mode and moves it by less than one grid cell; (b) the posterior median for every node, as the root already uses, which changes the meaning of the date; (c) no change, with the sensitivity stated in the documentation of the date outputs.

## Related

- [M-timetree-edge-message-regridding-limits-accuracy.md](M-timetree-edge-message-regridding-limits-accuracy.md)
