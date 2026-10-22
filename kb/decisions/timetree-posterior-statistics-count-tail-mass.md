# Posterior statistics of node times count the tail mass

**Type**: Intentional divergence from v0 and augur.

The root date, the survival functions of the coalescent prior, the 90% date intervals (`num_date_confidence`), the mass windows and the log marginal likelihood of node times all read one probability-mass model of a time distribution. The model counts the mass that the tail laws of a distribution carry beyond its grid. v0, and augur through v0's `get_max_posterior_region()`, integrate the grid samples only.

## Problem

A time distribution in v1 is a grid of negative log densities with a tail law on each side (`enum BoundaryBehavior`, [distribution-tails-and-arithmetic.md](distribution-tails-and-arithmetic.md)): a soft `Linear` side continues past the grid edge as an exponential, a `HardApproach` side follows a power law down to its hard bound. Before this decision, `fn log_mass()` counted the tail mass, while the posterior median, the survival functions and the 90% region integrated the grid alone. Example: density 1 on the grid $[-1, 0]$, a hard right side, and the left tail $e^{t + 1}$ for $t < -1$ hold the masses 1 on the grid and 1 in the tail, so the median is $-1$; the grid-only median was $-0.5$.

## Mass model

`struct MassProfile` with `struct MassFractions` and `enum TailMass` [packages/treetime-distribution/src/distribution_ops/mass_profile.rs](../../packages/treetime-distribution/src/distribution_ops/mass_profile.rs) is the only code that integrates a `DistributionFunction`. It returns the total mass, the mass between two points (also as a negative log, finite far below the peak), the CDF, the survival function, the quantile, the inverse of the survival function, and the density at a point. Masses are relative to the peak density, and the profile exposes the peak ordinate, so callers combine it with the amplitude in log space.

- **Regular cells**: the trapezoid rule of the grid samples. A cell of width $h$ between densities $p_0$ and $p_1$ holds $h\,(p_0 + p_1)/2$, and the density is linear in probability inside the cell, so partial masses and inverses are those of a linear density. This is the quadrature of v0 and of the earlier grid-only CDF. The convolution of two functions integrates each overlap by the trapezoid rule, and its two extreme outputs hold $h\,f_0 g_0/4$ and $h\,f_{F-1} g_{G-1}/4$, so the trapezoid total of its output equals the product of the trapezoid totals of its operands $f$ (samples $f_0, \dots, f_{F-1}$) and $g$ (samples $g_0, \dots, g_{G-1}$) up to $\tfrac{h^2}{4}\left(f_0 g_{G-1} + f_{F-1} g_0\right) - \tfrac{h^2}{8}\left(f_0 g_0 + f_{F-1} g_{G-1}\right)$
- **Hard cells**: a cell with one infinite ordinate (zero density) holds the mass of the power law that evaluation uses there (`fn GridFn::hard_cell()` [packages/treetime-grid/src/grid_fn.rs#L156](../../packages/treetime-grid/src/grid_fn.rs#L156)). With $p_a$ the density at the finite end and $\beta$ the exponent of the law, the cell holds $p_a\,h / (\beta + 1)$, and its partial mass at the distance $\delta$ from the hard end is $p_a\,h / (\beta + 1) \cdot (\delta / h)^{\beta + 1}$. A cell with two infinite ends holds no mass
- **Edge cells**: the outermost cell of a `HardEdgeCell` side holds the mass of its `EdgeCellLaw`, a power law in the distance to an origin beyond the edge times an exponential, integrated by Gauss-Legendre panels from its densest point ([kb/algo/distribution.md](../algo/distribution.md#edge-cell-with-an-offset-power-law)). Partial masses and the inverses of mass and density inside the cell follow the same law. A mode inside the cell is the peak of the posterior: the maximum posterior region grows from it, and its peak location is interior. No mass lies beyond the edge
- **Tails**: the closed-form masses of the `Linear` law (an exponential) and of the `HardApproach` power law, with their partial masses and inverses
- **Coordinates**: each cell takes the width between its end coordinates as the grid stores them after rounding, and positions inside a cell are offsets from these coordinates, so the model is continuous in the represented coordinate. On grids near the coordinate resolution this makes the mass depend on the rounding ([kb/issues/N-distribution-grids-near-the-coordinate-resolution-break-mass-and-survival-invariants.md](../issues/N-distribution-grids-near-the-coordinate-resolution-break-mass-and-survival-invariants.md)). At a hard point a density query returns the one-sided limit from inside the support
- **Numerics**: scaled and compensated sums, and a segment tree of cell masses for intervals far below the peak and across valleys, so small masses keep their digits; each inverse starts from the end at which its target accumulates
- **Undefined total**: an `Error` side or a soft tail that does not decay has no finite mass. `TailMass::Required` then reports no fractions; `TailMass::ZeroWhereUndefined` counts no mass beyond such a side

### Rejected quadrature

The exact integral of the log-linear interpolant in each cell undercounts a Gaussian by about $h^2 / (12\,\sigma^2)$ relative, where $\sigma$ is its standard deviation: $4.6 \cdot 10^{-5}$ on 300 points over $\pm 3.5\,\sigma$, against $2.8 \cdot 10^{-7}$ for the trapezoid. It also breaks the identity that an FFT convolution adds the log masses of its operands. A zero-mass hard cell with a log-linear limit interpolation was rejected too: a query one rounding error past a stored point next to such a cell read an infinite neg-log and left 46 nodes without a date.

## Statistics on the model

`packages/treetime/src/timetree/max_posterior_region.rs`:

- **Root date**: `fn posterior_median()` takes the quantile 0.5 under `TailMass::ZeroWhereUndefined` ([timetree-root-time-at-posterior-median.md](timetree-root-time-at-posterior-median.md))
- **90% interval**: `fn max_posterior_region()` finds the density level whose region holds 90% of the model mass. Each side of the peak maps a level to its first crossing outward from the peak, through the running maximum of the side ([kb/v0-errata/timetree-max-posterior-region-multimodal-side.md](../v0-errata/timetree-max-posterior-region-multimodal-side.md)). Inside a regular cell the crossing is linear in probability, inside a hard cell it inverts the power law, and in a tail it inverts the tail law (`fn SoftTailLaw::position_at()`, `fn HardApproachLaw::position_at()`). Brent's method (`argmin`) solves for the level only inside one continuous level interval; at a plateau, a zero-density gap or a flat tail the mass jumps, and the region ends inside the plateau where it holds 90%. The solve runs in coordinates centered on the peak, and every region is checked to hold 90% of the mass within rounding. A peak at a grid end whose side carries no tail mass gives `[Q(0.1), t_max]` or `[t_min, Q(0.9)]`, where $Q$ is the quantile
- **Survival functions**: `fn posterior_survival()` represents the survival function $S(t) = P(T > t)$ of a node posterior by 55 knots at fixed probability levels (`fn survival_levels()` in [packages/treetime/src/coalescent/mean_field.rs](../../packages/treetime/src/coalescent/mean_field.rs)): $10^{-9}$, $10^{-7}$, $10^{-5}$, $10^{-4}$, $10^{-3}$, $3 \cdot 10^{-3}$, $10^{-2}$, $2 \cdot 10^{-2}$ and $3.5 \cdot 10^{-2}$ in each tail and steps of 0.025 between. Each knot lies at the model quantile of its level (`fn MassFractions::quantile()`, `fn MassFractions::survival_inverse()`) and stores the model survival there, and $S$ is linear between knots. [Error of the survival levels](#error-of-the-survival-levels) derives the interpolation error. The Gaussian seed of the init phase uses the same levels on its truncated normal. The extreme levels must give finite times; a quantile that is not finite is an error that names the tail slope. At the resolution of the calendar coordinates a rounded cutoff can omit more mass ([kb/issues/N-distribution-grids-near-the-coordinate-resolution-break-mass-and-survival-invariants.md](../issues/N-distribution-grids-near-the-coordinate-resolution-break-mass-and-survival-invariants.md))

  > [!WARNING]
  > **Needs review.** The validation numbers of the next paragraph carry no measurement provenance (commands, data sets, code states) and have not been reproduced: [kb/issues/M-timetree-low-signal-rate-bias-claims-not-reproduced.md](../issues/M-timetree-low-signal-rate-bias-claims-not-reproduced.md).

  A cubic Hermite interpolant with an interpolation error below $10^{-8}$ needs a median of about 290 knots per posterior. The coalescent prior does not need that precision: compared with such an interpolant on the validation grid (240 sets, `--coalescent-opt` and `--coalescent 10`, 6 rounds) the linear survival on fixed levels gives the same clock rates, dates and root dates to three digits wherever there is signal, the same 90% interval coverage (0.898) on the calibration runs, and 8% shorter run times with 17% less memory; the full timetree suite (2970 runs) changes median date RMSE and rate error by at most $10^{-3}$
- **Ranges and points**: a `Range` posterior keeps closed forms (median at the midpoint, linear survival, region at the latest 90% of the range); a `Point` posterior gives its time and the region `[t, t]`

The coalescent prior integrates the expected lineage counts that these survival functions give, with $1/T_c$ constant on each cell of its grid, by the trapezoid rule over the linear pieces, split at the cells (`struct CoalescentPrior`, [packages/treetime/src/coalescent/mean_field.rs](../../packages/treetime/src/coalescent/mean_field.rs)).

## Error of the survival levels

The survival levels replace the exact model survival $S(t)$ by its linear interpolant $L(t)$ through the knots $(t_k, S(t_k))$, with $L = 1$ before the first knot and $L = 0$ from the last knot on. Symbols:

- $c_k = 1 - S(t_k)$: the level of knot $k$ as a CDF value
- $\Delta t_k = t_{k+1} - t_k$: the width of the piece between knots $k$ and $k + 1$
- $p(t) = -S'(t)$: the normalized posterior density, $p'$ its slope
- $r = c_{k+1} / c_k$: the ratio of two neighboring levels

Three bounds hold on each piece:

- **Bracket**: $S$ and $L$ both lie between $S(t_k)$ and $S(t_{k+1})$, so $\lvert S - L \rvert \le c_{k+1} - c_k$: 0.025 between interior levels, at most 0.015 between the tail levels from $10^{-2}$ on, and at most $10^{-5}$ below the level $10^{-5}$. This bound also holds across a gap of zero density, where $S$ is flat and $L$ is not
- **Smooth density**: $\lvert S - L \rvert \le \tfrac18 \Delta t_k^2 \max \lvert p' \rvert$, the error of linear interpolation of a function with second derivative $-p'$
- **Exponential tail**: on a `Linear` soft tail the CDF is exponential in $t$, and the largest interpolation error of the piece is $c_k\,(1 + u (r - 1) - r^u)$ with $r^u = (r - 1)/\ln r$, independent of the tail rate. The largest value over the tail levels is $1.26 \cdot 10^{-3}$, on the piece between the levels 0.05 and 0.075 ($r = 1.5$); the pieces $2 \cdot 10^{-2}$ to $3.5 \cdot 10^{-2}$ and $3 \cdot 10^{-3}$ to $10^{-2}$ give $1.04 \cdot 10^{-3}$ and $1.03 \cdot 10^{-3}$

The own integral of a node, $\int L(t)/T_c(t)\,dt$, is exact for the linear pieces. Against $\int S/T_c$, each piece on which $p'$ is constant contributes $-\Delta t^3 p' / (12\,T_c)$, so pieces with rising and with falling density cancel in part.

### Bounds at the asserted points of the tests

The tests compare the representation with closed forms. An independent script builds the knots from the closed form of each fixture (knot at the exact quantile of its level, value of the exact survival there) and computes $\lvert S - L \rvert$ at the asserted points. The code gives the same errors to at least 4 digits, so no test sees an error beyond the representation. Each literal is the tightest passing `1e-N` above the derived error.

Tests of `packages/treetime/src/coalescent/__tests__/test_mean_field.rs`:

- `test_mean_field_segment_statistics_read_the_survival_inside_a_long_earlier_tail`: $S(-2)$ of a unit box on $[-1, 0]$ with the tail $e^{t+1}$ (CDF $e^{t+1}/2$ in the tail). The point $-2$ lies between the knots at the CDF levels 0.175 and 0.2, $\Delta t = \ln(8/7) = 0.1335$. Error $3.88 \cdot 10^{-4}$; smooth-density bound $\tfrac18 \Delta t^2 \cdot 0.2 = 4.46 \cdot 10^{-4}$. Literal `1e-3` in both assertions
- `test_mean_field_segment_statistics_read_the_survival_of_a_narrow_posterior_with_a_long_soft_tail`: a tent density on $[0, 1]$ (from 0.1 to 1 and back) with a soft tail of rate 0.2 to the left that holds 48% of the mass, read at 2220 boundaries from $-110$ to 1. The largest error, $3.69 \cdot 10^{-3}$ at $t = 0.05$, lies on the piece $[-0.0125, 0.1201]$ between the levels 0.475 and 0.5, which spans the junction of the tail and the steep rising flank of the tent ($p' = 1.8/1.05$); smooth-density bound $3.77 \cdot 10^{-3}$. The tail pieces stay below $1.26 \cdot 10^{-3}$. Literal `1e-2`
- `test_mean_field_own_integral_follows_a_tc_step_inside_a_knot_interval`: own integral of a piecewise-linear density (samples 1, 3, 2, 4, 1 on $[0, 4]$) across a $T_c$ step from 2 to 0.5 at $t = 0.5$. Errors $1.4 \cdot 10^{-5}$, $1.42 \cdot 10^{-4}$, $1.46 \cdot 10^{-4}$, $2.15 \cdot 10^{-4}$ and $5.3 \cdot 10^{-5}$ at $x = 0.25, 0.75, 1.5, 2.5, 3.99$. The sum of the absolute piece errors reaches $1.0 \cdot 10^{-3}$ at $x = 3.99$; the signed sum is 19 times smaller, because the density rises and falls in turn. Literal `1e-3`
- `test_mean_field_gaussian_seed_truncates_at_the_earliest_included_tip`: the truncated normal (mean 2005, standard deviation 1, upper bound 2010) read at 2000, between the knots at the CDF levels $10^{-7}$ ($t = 1999.80$) and $10^{-5}$ ($t = 2000.73$). $1 - L(2000) = 2.21 \cdot 10^{-6}$, of which $2.87 \cdot 10^{-7}$ is the mass of the normal below 2000; bracket bound $10^{-5}$. Literal `1e-5`

Tests of `packages/treetime/src/timetree/__tests__/test_posterior_statistics.rs`, largest error over 17 points per piece:

- Unit box with an earlier tail of mass 1, and a narrow box with a later tail of rate 1: $1.26 \cdot 10^{-3}$, the exponential-tail maximum. Literal `1e-2`
- Box with a steep later tail of rate 100: $1.03 \cdot 10^{-3}$. Literal `1e-2`
- Box with a `HardApproach` tail of exponent 1: $6.31 \cdot 10^{-4}$; with exponent 2.5: $9.02 \cdot 10^{-4}$. Literal `1e-3`
- Two boxes around a gap of zero density: $2.19 \cdot 10^{-2}$, the piece across the gap; bracket bound 0.025. Literal `1e-1`
- Cells $10^{-10}$ wide at calendar time 2000: $4.6 \cdot 10^{-4}$ with trapezoid cells, rounding with a flat cell. Literals `1e-3` and `1e-15`
- A uniform density whose survival is linear: the representation is exact, the error is rounding ($1.4 \cdot 10^{-17}$). Literal `1e-15`
- The knot times and the stored survival values match the closed forms to rounding: the knot-time literals are those of the cubic Hermite tests, and the knot values use `1e-15`

### Effect on the exact Kingman posterior

`test_marginal_exact_mean_field_prior_matches_the_kingman_posterior` ([packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs)) compares node CDFs under the mean-field prior with the exact Kingman posterior summed on a 0.02-year grid. Halving the step of that reference changes its CDFs by $1.1 \cdot 10^{-5}$ to $7.1 \cdot 10^{-4}$ per case and node. Against a cubic Hermite survival function with an error below $10^{-8}$ and the same convolution, the survival levels change the largest CDF difference per node by $3 \cdot 10^{-7}$ to $8.5 \cdot 10^{-5}$. The change exceeds the reference change at one node only, R of the overlapping clades: from $8.566 \cdot 10^{-3}$ to $8.515 \cdot 10^{-3}$ with $T_c = 1$ (reference change $2.1 \cdot 10^{-5}$), and from $1.512 \cdot 10^{-3}$ to $1.524 \cdot 10^{-3}$ with $T_c = 10$ (reference change $1.1 \cdot 10^{-5}$). Both changes are below 1% of the difference to the exact posterior, which the mean-field approximation sets.

## Clip of the interval

`fn node_confidence_interval()` clips an interval only at the hard sides of the node's subtree distribution (the grid edge of a `Hard` or `Error` side, the bound of a `HardApproach` side, no clip on a soft side), never at the grid end of the posterior. Why the subtree distribution and not the posterior: [timetree-forward-posterior-trimming.md](timetree-forward-posterior-trimming.md#confidence-intervals-are-clipped-at-physical-bounds).

## Effect

Measured with `timetree --time-marginal always` against the grid-only statistics, before the forward-posterior trim: root dates change by at most $1.5 \cdot 10^{-7}$ years and 90% bounds by at most $7 \cdot 10^{-5}$ years on `data/flu/h3n2/20`, `data/ebola/20`, and `data/flu/h3n2/20` with `--coalescent opt`. On these runs the tails hold little mass. The difference grows where a posterior holds noticeable mass in a tail law, for example after a mass window cut near the peak or at a soft side close to the peak.

The region now holds 90% of the mass on every posterior. Before, a posterior whose cells next to a zero-density point held more than 10% of its mass fell back to the bounds of its finite support, with less than 90% inside (grid $t = [0, 1, 2, 3, 4, 5]$, neg-logs $[1, 0, 5, 0.5, 0.5, \infty]$: $[0, 4]$ holding 0.874).

## Divergence from v0 and augur

v0 builds its CDF by the trapezoid rule of the grid samples (`_ml_t_marginal`, [packages/legacy/treetime/treetime/clock_tree.py#L931-L949](../../packages/legacy/treetime/treetime/clock_tree.py#L931-L949)) and searches the 90% level with Brent on $(0.9 - \text{mass}(h))^2$ ([packages/legacy/treetime/treetime/clock_tree.py#L1146-L1230](../../packages/legacy/treetime/treetime/clock_tree.py#L1146-L1230)), where $h$ is the level. augur writes `num_date_confidence` from it. v1 and v0 therefore agree on regular cells and differ by the mass of tails and hard cells, and where v0's Brent search stops on a plateau ([kb/v0-errata/timetree-max-posterior-region-brent-plateau.md](../v0-errata/timetree-max-posterior-region-brent-plateau.md)).

## Code references

- `struct MassProfile`, `struct MassFractions`, `enum TailMass` [packages/treetime-distribution/src/distribution_ops/mass_profile.rs](../../packages/treetime-distribution/src/distribution_ops/mass_profile.rs)
- `fn log_mass()`, `fn mass_bounded_domain()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs)
- `fn GridFn::hard_cell()`, `fn GridFn::hard_cell_fit_range()` [packages/treetime-grid/src/grid_fn.rs#L156-L190](../../packages/treetime-grid/src/grid_fn.rs#L156-L190)
- `fn max_posterior_region()`, `fn posterior_median()`, `fn posterior_survival()` [packages/treetime/src/timetree/max_posterior_region.rs](../../packages/treetime/src/timetree/max_posterior_region.rs)
- `fn NodePosterior::summarize()` [packages/treetime/src/timetree/inference/result.rs#L109](../../packages/treetime/src/timetree/inference/result.rs#L109)
