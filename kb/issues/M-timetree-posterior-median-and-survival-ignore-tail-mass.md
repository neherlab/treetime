# Posterior median and survival function ignore the tail mass of a distribution

## Summary

`fn posterior_median()` and `fn posterior_survival()` build their cumulative distribution (`struct GridCdf`) from the stored grid of a time distribution only. A function distribution can also carry probability mass outside its grid: its left and right boundary behaviors (`BoundaryBehavior::Linear` soft tails, `BoundaryBehavior::HardApproach`) describe how the density continues past the first and the last grid point. `fn log_mass()` counts this tail mass; the median and the survival function do not. When a tail holds a noticeable part of the mass, both results are computed from the wrong distribution.

Locations:

- `pub(crate) fn posterior_median()` [packages/treetime/src/timetree/max_posterior_region.rs#L29-L41](../../packages/treetime/src/timetree/max_posterior_region.rs#L29-L41)
- `pub(crate) fn posterior_survival()` [packages/treetime/src/timetree/max_posterior_region.rs#L43-L57](../../packages/treetime/src/timetree/max_posterior_region.rs#L43-L57)
- `struct GridCdf` (trapezoid integration over the grid points) [packages/treetime/src/timetree/max_posterior_region.rs#L208-L232](../../packages/treetime/src/timetree/max_posterior_region.rs#L208-L232)
- Tail-aware mass: `pub(crate) fn mass_profile()` adds `left_mass` and `right_mass` from `fn tail_mass()` to the grid mass [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L226-L283](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L226-L283), used by `pub(crate) fn log_mass()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L10-L21](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L10-L21)

## Example

Take a distribution with density $p(t) = 1$ on the grid $[-1, 0]$, a hard right boundary, and the soft left tail $p(t) = e^{t + 1}$ for $t < -1$. The grid holds mass $1$ and the tail holds $\int_{-\infty}^{-1} e^{t+1}\,dt = 1$, so the total mass is $2$ and the median is the point with cumulative mass $1$, which is $t = -1$. `fn posterior_median()` integrates the grid only and returns $-0.5$.

## Consumers

- **Root date.** The forward pass commits the median of the root posterior as the root time: `fn propagate_distributions_forward_node()` calls `posterior_median` [packages/treetime/src/timetree/inference/forward_pass.rs#L94-L101](../../packages/treetime/src/timetree/inference/forward_pass.rs#L94-L101)
- **Coalescent prior.** `fn LineageExpectation::from_posterior()` builds the survival function $S_u(t) = P(t_u > t)$ of each node from `posterior_survival` [packages/treetime/src/coalescent/mean_field.rs#L165-L186](../../packages/treetime/src/coalescent/mean_field.rs#L165-L186). The expected lineage counts $\bar k(t) = \sum_u w_u S_u(t)$, the merger rates of polytomy resolution, the $T_c$ fit, and the mean-field prior all read these survival functions
- **Related.** `pub fn max_posterior_region()` also builds its region from the grid only [packages/treetime/src/timetree/max_posterior_region.rs#L10-L27](../../packages/treetime/src/timetree/max_posterior_region.rs#L10-L27), so a confidence interval leaves out the tail mass in the same way

> [!IMPORTANT]
> **Investigation required.** The size of the effect on real data is unmeasured. Collect, for the root posterior and the internal-node posteriors of representative runs (for example `data/flu/h3n2/500`, `data/ebola/362`, `data/mpox/clade-ii/500` with `timetree --coalescent-opt`), the fraction of the mass that lies in the tails, and the shift of the root median and of $\bar k(t)$ when the tails are counted.

## Fix direction

Give the distribution crate one tail-aware cumulative mass, next to `fn mass_profile()` in the mass model, that integrates the grid and adds the tail masses of `fn tail_mass()` at the ends. Use it in `fn log_mass()`, in the posterior median, in the survival function, and in the maximum posterior region, so all four read the same mass. A quantile that falls in a tail needs the inverse of the tail law (`SoftTailLaw`, `HardApproachLaw`), and a survival function on the grid starts below 1 by the left tail mass and ends above 0 by the right tail mass.
