# Grids near the coordinate resolution break mass and survival invariants

Time distributions live on calendar coordinates, where the unit in the last place (ulp) of a date near 2000 is $u = 2^{-42} \approx 2.3 \cdot 10^{-13}$ years. When a grid spacing or a tail decay scale comes close to this resolution, the mass model (`struct MassProfile` in `packages/treetime-distribution/src/distribution_ops/mass_profile.rs`) and the survival functions of the coalescent prior lose three invariants. Each cell of the mass model has the width between its end coordinates as the grid stores them after rounding ([packages/treetime-distribution/src/distribution_ops/mass_profile.rs#L62](../../packages/treetime-distribution/src/distribution_ops/mass_profile.rs#L62)). This rule causes two of the three defects below.

## Zero-width cells

When the grid spacing is below half an ulp, neighboring grid points round to the same coordinate and their cell has zero width and zero mass. A grid whose points all fall on one coordinate has no mass at all, and the survival functions (`fn posterior_survival()` in `packages/treetime/src/timetree/max_posterior_region.rs`) then stop with an error.

## Survival cutoff rounded to the grid edge

`fn function_survival()` takes the $10^{-9}$ quantiles as the ends of the survival interpolant and does not check the mass that remains beyond them after rounding ([packages/treetime/src/timetree/max_posterior_region.rs#L154](../../packages/treetime/src/timetree/max_posterior_region.rs#L154)); `fn NodeSurvival::at()` returns 0 at and after the last knot ([packages/treetime/src/coalescent/mean_field.rs#L646](../../packages/treetime/src/coalescent/mean_field.rs#L646)). For a tail whose decay scale is below the coordinate resolution, the quantile rounds to the grid edge. Case: start 2000, spacing $10^{-10}$, neg-logs $[0, 0]$, a `Hard` left side and a `Linear` right side of slope $10^{15}$:

- the stored end point lies $H = 1.000444 \cdot 10^{-10}$ after the start
- the exact survival there is $1 / (1 + 10^{15} H) = 9.995 \cdot 10^{-6}$
- the upper quantile lies $9.2 \cdot 10^{-15}$ past the end point, which rounds to the end point, and the interpolant reads 0 there

The error is about $10^{-5}$. Since 2026-10-21 the survival functions interpolate linearly between fixed probability levels, with an error bound of the step between levels ($\ge 2 \cdot 10^{-9}$ in the far tails), so the error here no longer exceeds the stated bound by orders of magnitude, but the cutoff still omits the mass beyond the rounded end.

## Mass changes under negation and translation

`fn GridFn::negate_arg_inplace()` rebuilds a grid from its negated end point and its stored spacing ([packages/treetime-grid/src/grid_fn.rs#L375](../../packages/treetime-grid/src/grid_fn.rs#L375)), and the rounded coordinates of the new grid can differ from the reflected old ones. Case: start 2000, spacing $1.5\,u$, neg-logs $[0, 0, \ln 4, \ln 4]$, both sides `Hard`:

| Quantity | Original | After negation |
| --- | --- | --- |
| Coordinates relative to the start, in units of $u$ | $[0, 2, 3, 4]$ | $[-4, -2, -1, 0]$ |
| Densities | $[1, 1, 0.25, 0.25]$ | $[0.25, 0.25, 1, 1]$ |
| Total mass | $2.875\,u$ | $2.125\,u$ |

The log mass changes by 0.302. The confidence region translates the posterior to coordinates centered on its peak ([packages/treetime/src/timetree/max_posterior_region.rs#L359](../../packages/treetime/src/timetree/max_posterior_region.rs#L359)), so its mass depends on the coordinates in the same way. A mass model that integrates widths of the nominal spacing gives both grids $2.8125\,u$, but then the survival interpolant jumps at the same scale: exact survival and exact mass conservation under these transformations cannot both hold on rounded coordinates.

> [!IMPORTANT]
> **Investigation required.** No input is known to produce such grids. Log the smallest grid spacing and the steepest tail slope, relative to the coordinate ulp, of every node posterior and message on the example datasets, the validation subset, and runs at fixed clock rates far above the estimated ones (for example `data/mpox/clade-ii/1000` at `--clock-rate=0.5`).

## Fix direction

One owner policy for distributions near the coordinate resolution: either represent a distribution narrower than the resolution as a point at its coordinate, with its total mass as the amplitude, or compute masses, quantiles, negation and translation in local coordinates relative to one origin throughout. Check the mass that a rounded cutoff omits, and return an error when the coordinate resolution prevents the survival bound.

## Validation

- A distribution narrower than one ulp gives the median, interval and survival function of a point, without an error
- The survival case above stays within $10^{-8}$ through `fn NodeSurvival::of_posterior()`, on both reflected sides
- Negation and translation of the case above keep the total mass within rounding
