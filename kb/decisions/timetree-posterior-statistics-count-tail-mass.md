# Posterior statistics of node times count the tail mass

**Type**: Intentional divergence from v0 and augur.

The root date, the survival functions of the coalescent prior, the 90% date intervals (`num_date_confidence`), the mass windows and the log marginal likelihood of node times all read one probability-mass model of a time distribution. The model counts the mass that the tail laws of a distribution carry beyond its grid. v0, and augur through v0's `get_max_posterior_region()`, integrate the grid samples only.

## Problem

A time distribution in v1 is a grid of negative log densities with a tail law on each side (`enum BoundaryBehavior`, [distribution-tails-and-arithmetic.md](distribution-tails-and-arithmetic.md)): a soft `Linear` side continues past the grid edge as an exponential, a `HardApproach` side follows a power law down to its hard bound. Before this decision, `fn log_mass()` counted the tail mass, while the posterior median, the survival functions and the 90% region integrated the grid alone. Example: density 1 on the grid $[-1, 0]$, a hard right side, and the left tail $e^{t + 1}$ for $t < -1$ hold the masses 1 on the grid and 1 in the tail, so the median is $-1$; the grid-only median was $-0.5$.

## Mass model

`struct MassProfile` with `struct MassFractions` and `enum TailMass` [packages/treetime-distribution/src/distribution_ops/mass_profile.rs](../../packages/treetime-distribution/src/distribution_ops/mass_profile.rs) is the only code that integrates a `DistributionFunction`. It returns the total mass, the mass between two points (also as a negative log, finite far below the peak), the CDF, the survival function, the quantile, the inverse of the survival function, and the density at a point. Masses are relative to the peak density, and the profile exposes the peak ordinate, so callers combine it with the amplitude in log space.

- **Regular cells**: the trapezoid rule of the grid samples. A cell of width $h$ between densities $p_0$ and $p_1$ holds $h\,(p_0 + p_1)/2$, and the density is linear in probability inside the cell, so partial masses and inverses are those of a linear density. This is the quadrature of v0 and of the earlier grid-only CDF, and its total is the one that the Riemann sums of the FFT convolution conserve
- **Hard cells**: a cell with one infinite ordinate (zero density) holds the mass of the power law that evaluation uses there (`fn GridFn::hard_cell()` [packages/treetime-grid/src/grid_fn.rs#L156](../../packages/treetime-grid/src/grid_fn.rs#L156)). With $p_a$ the density at the finite end and $\beta$ the exponent of the law, the cell holds $p_a\,h / (\beta + 1)$, and its partial mass at the distance $\delta$ from the hard end is $p_a\,h / (\beta + 1) \cdot (\delta / h)^{\beta + 1}$. A cell with two infinite ends holds no mass
- **Tails**: the closed-form masses of the `Linear` law (an exponential) and of the `HardApproach` power law, with their partial masses and inverses
- **Numerics**: scaled and compensated sums, and a segment tree of cell masses for intervals far below the peak and across valleys, so small masses keep their digits; each inverse starts from the end at which its target accumulates
- **Undefined total**: an `Error` side or a soft tail that does not decay has no finite mass. `TailMass::Required` then reports no fractions; `TailMass::ZeroWhereUndefined` counts no mass beyond such a side

### Rejected quadrature

The exact integral of the log-linear interpolant in each cell undercounts a Gaussian by about $h^2 / (12\,\sigma^2)$ relative, where $\sigma$ is its standard deviation: $4.6 \cdot 10^{-5}$ on 300 points over $\pm 3.5\,\sigma$, against $2.8 \cdot 10^{-7}$ for the trapezoid. It also breaks the identity that an FFT convolution adds the log masses of its operands. A zero-mass hard cell with a log-linear limit interpolation was rejected too: a query one rounding error past a stored point next to such a cell read an infinite neg-log and left 46 nodes without a date.

## Statistics on the model

`packages/treetime/src/timetree/max_posterior_region.rs`:

- **Root date**: `fn posterior_median()` takes the quantile 0.5 under `TailMass::ZeroWhereUndefined` ([timetree-root-time-at-posterior-median.md](timetree-root-time-at-posterior-median.md))
- **90% interval**: `fn max_posterior_region()` finds the density level whose region holds 90% of the model mass. Each side of the peak maps a level to its first crossing outward from the peak, through the running maximum of the side ([kb/v0-errata/timetree-max-posterior-region-multimodal-side.md](../v0-errata/timetree-max-posterior-region-multimodal-side.md)). Inside a regular cell the crossing is linear in probability, inside a hard cell it inverts the power law, and in a tail it inverts the tail law (`fn SoftTailLaw::position_at()`, `fn HardApproachLaw::position_at()`). Brent's method (`argmin`) solves for the level only inside one continuous level interval; at a plateau, a zero-density gap or a flat tail the mass jumps, and the region ends inside the plateau where it holds 90%. The solve runs in coordinates centered on the peak, and every region is checked to hold 90% of the mass within rounding. A peak at a grid end whose side carries no tail mass gives `[Q(0.1), t_max]` or `[t_min, Q(0.9)]`, where $Q$ is the quantile
- **Survival functions**: `fn posterior_survival()` interpolates the survival function $S(t) = P(T > t)$ of a node posterior by cubic Hermite pieces whose derivatives are the exact model density (`fn MassFractions::density()`). Knots lie on the grid points, where the interpolation is exact in regular cells because $S$ is quadratic there; hard cells and tails get extra knots until an a-priori bound on the interpolation error, from the fourth derivative of the cell or tail law, holds. The pieces run from the $10^{-9}$ quantile to the upper $10^{-9}$ quantile (`SURVIVAL_TRIM`), so they continue into the tails; the error bound is $10^{-8}$ (`SURVIVAL_INTERPOLATION_ERROR`, [packages/treetime/src/coalescent/mean_field.rs#L21-L22](../../packages/treetime/src/coalescent/mean_field.rs#L21-L22)). A posterior needs a median of about 290 knots. A quantile that is not finite is an error that names the tail slope
- **Ranges and points**: a `Range` posterior keeps closed forms (median at the midpoint, linear survival, region at the latest 90% of the range); a `Point` posterior gives its time and the region `[t, t]`

The coalescent prior integrates the expected lineage counts that these survival functions give, with $1/T_c$ constant on each cell of its grid, in closed form over the Hermite pieces (`struct CoalescentPrior`, [packages/treetime/src/coalescent/mean_field.rs](../../packages/treetime/src/coalescent/mean_field.rs)).

## Clip of the interval

`fn node_confidence_interval()` clips an interval only at the hard sides of the node's subtree distribution (the grid edge of a `Hard` or `Error` side, the bound of a `HardApproach` side, no clip on a soft side), never at the grid end of the posterior. Why the subtree distribution and not the posterior: [timetree-forward-posterior-trimming.md](timetree-forward-posterior-trimming.md#confidence-intervals-are-clipped-at-physical-bounds).

## Effect

Measured with `timetree --time-marginal always` against the grid-only statistics, before the forward-posterior trim: root dates change by at most $1.5 \cdot 10^{-7}$ years and 90% bounds by at most $7 \cdot 10^{-5}$ years on `data/flu/h3n2/20`, `data/ebola/20`, and `data/flu/h3n2/20` with `--coalescent-opt`. On these runs the tails hold little mass. The difference grows where a posterior holds noticeable mass in a tail law, for example after a mass window cut near the peak or at a soft side close to the peak.

The region now holds 90% of the mass on every posterior. Before, a posterior whose cells next to a zero-density point held more than 10% of its mass fell back to the bounds of its finite support, with less than 90% inside (grid $t = [0, 1, 2, 3, 4, 5]$, neg-logs $[1, 0, 5, 0.5, 0.5, \infty]$: $[0, 4]$ holding 0.874).

## Divergence from v0 and augur

v0 builds its CDF by the trapezoid rule of the grid samples (`_ml_t_marginal`, [packages/legacy/treetime/treetime/clock_tree.py#L931-L949](../../packages/legacy/treetime/treetime/clock_tree.py#L931-L949)) and searches the 90% level with Brent on $(0.9 - \text{mass}(h))^2$ ([packages/legacy/treetime/treetime/clock_tree.py#L1146-L1230](../../packages/legacy/treetime/treetime/clock_tree.py#L1146-L1230)), where $h$ is the level. augur writes `num_date_confidence` from it. v1 and v0 therefore agree on regular cells and differ by the mass of tails and hard cells, and where v0's Brent search stops on a plateau ([kb/v0-errata/timetree-max-posterior-region-brent-plateau.md](../v0-errata/timetree-max-posterior-region-brent-plateau.md)).

## Code references

- `struct MassProfile`, `struct MassFractions`, `enum TailMass` [packages/treetime-distribution/src/distribution_ops/mass_profile.rs](../../packages/treetime-distribution/src/distribution_ops/mass_profile.rs)
- `fn log_mass()`, `fn mass_bounded_domain()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs)
- `fn GridFn::hard_cell()`, `fn GridFn::hard_cell_fit_range()` [packages/treetime-grid/src/grid_fn.rs#L156-L190](../../packages/treetime-grid/src/grid_fn.rs#L156-L190)
- `fn max_posterior_region()`, `fn posterior_median()`, `fn posterior_survival()` [packages/treetime/src/timetree/max_posterior_region.rs](../../packages/treetime/src/timetree/max_posterior_region.rs)
- `fn NodePosterior::summarize()` [packages/treetime/src/timetree/inference/result.rs#L109](../../packages/treetime/src/timetree/inference/result.rs#L109)
