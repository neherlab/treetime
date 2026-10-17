# Maximum posterior region holds less than the requested mass when no finite level reaches it

## Summary

`fn highest_density_interval()` [packages/treetime/src/timetree/max_posterior_region.rs](../../packages/treetime/src/timetree/max_posterior_region.rs) raises a level above the peak until the region between the two level crossings holds the requested fraction of the posterior mass (90% for `timetree --confidence`). When no finite level reaches the fraction, it returns the bounds of `fn LevelSide::support_end()`: on each side of the peak, the last grid point before the first point with infinite negative log-likelihood (zero density). This region can hold less than the requested fraction.

The cause is that two parts of the function discretize the density differently between a finite point and an infinite point:

- The level map (`struct LevelSide`) interpolates the negative log-likelihood linearly, so a cell that ends at an infinite value has zero density on its whole open interval. No finite level reaches into it
- The CDF (`struct GridCdf`) integrates the density with the trapezoid rule, so the same cell holds `(t1 - t0) * p0 / 2` of mass

The mass of these end cells is outside every region that the level search can return.

## Reproduction

Grid `t = [0, 1, 2, 3, 4, 5]`, negative log-likelihoods `[1, 0, 5, 0.5, 0.5, inf]`, fraction 0.9 (test `test_max_posterior_region_spans_support_when_no_level_holds_fraction` in [packages/treetime/src/timetree/**tests**/test_max_posterior_region.rs](../../packages/treetime/src/timetree/__tests__/test_max_posterior_region.rs)):

- Trapezoid cell masses: $\frac{e^{-1}+1}{2}$, $\frac{1+e^{-5}}{2}$, $\frac{e^{-5}+e^{-0.5}}{2}$, $e^{-0.5}$, $\frac{e^{-0.5}}{2}$; total $\approx 2.4037$
- v1 returns `[0, 4]`. It holds $1 - \frac{e^{-0.5}/2}{2.4037} \approx 0.874$ of the mass, below the requested 0.9. The test asserts these bounds without asserting their mass
- v0 `ClockTree.get_max_posterior_region()` ([packages/legacy/treetime/treetime/clock_tree.py#L1146-L1230](../../packages/legacy/treetime/treetime/clock_tree.py#L1146-L1230)), run directly on the same posterior with the CDF that `_ml_t_marginal` builds ([packages/legacy/treetime/treetime/clock_tree.py#L931-L949](../../packages/legacy/treetime/treetime/clock_tree.py#L931-L949)), returns `[0, 2]`, which holds 0.494 of the mass by v0's own CDF. Its Brent search ends on the plateau above the highest finite height (objective 0.165, reported as success). The second mode at `t = 3..4` is lost through the non-monotone side described in [kb/v0-errata/timetree-max-posterior-region-multimodal-side.md](../v0-errata/timetree-max-posterior-region-multimodal-side.md)

v0 gives no target value for this case: both v0 and v1 return less than the requested mass, and v0 returns less than v1. v0 and v1 use the same trapezoid CDF, so the mass of the end cells is the same in both.

> [!IMPORTANT]
> **Decision required.** v1 diverges from the docstring contract of v0 ("contains the specified fraction of the probability mass", [packages/legacy/treetime/treetime/clock_tree.py#L1149-L1150](../../packages/legacy/treetime/treetime/clock_tree.py#L1149-L1150)) and v0 does not meet it either, so the fix needs a decision on the target behavior. Options:
>
> - Extend the fallback bounds into the end cells by the CDF until the region holds the fraction. The region then holds exactly the requested mass; the bound inside an end cell is not a level crossing, because the level map gives that cell zero density
> - Return the whole grid, including the infinite end points, as the fallback. The region then holds all of the mass and overstates the uncertainty, as v0 does for shallow posteriors ([kb/v0-errata/timetree-max-posterior-region-brent-plateau.md](../v0-errata/timetree-max-posterior-region-brent-plateau.md))
> - Make the CDF agree with the level map (zero mass in a cell that ends at an infinite value). This changes the CDF on every posterior with an infinite end point, including the golden-master comparison with v0 `marginal_cdf`
> - Keep the current behavior and document that the region can hold less than the fraction on such posteriors

> [!IMPORTANT]
> **Investigation required.** How often the fallback runs on real data is unmeasured. It needs a posterior whose cells next to infinite points hold more than `1 - fraction` of the mass, which a fine time grid makes unlikely. Count the fallback cases per node in `timetree --confidence` runs on the example datasets.

## Impact

- `num_date_confidence` and the confidence interval outputs of `timetree --confidence` can be narrower than a 90% region for posteriors where the fallback runs
