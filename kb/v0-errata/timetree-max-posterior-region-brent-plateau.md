# get_max_posterior_region returns the whole support for shallow posteriors

## v0 location

`ClockTree.get_max_posterior_region()` (`#get_max_posterior_region`) [packages/legacy/treetime/treetime/clock_tree.py#L1213-L1226](../../packages/legacy/treetime/treetime/clock_tree.py#L1213-L1226)

## Erratum

For a posterior with its peak inside the grid, v0 searches for the level `h` above the peak at which the region holds 90% of the mass. It minimizes `(0.9 - mass(h))**2` with `scipy.optimize.minimize_scalar(..., bracket=[0, 10], method='brent')`. When both sides of the peak stay below height 10, every level above the higher side puts both boundaries on the ends of the grid, so `mass(h) = 1` and the objective is a flat `0.01`. The bracket search starts at `h = 10` on this plateau, moves further away from the peak, and Brent stops there. v0 then reports the whole support of the posterior. A level with mass exactly 0.9 exists below the plateau, because the mass is continuous in `h`.

Shallow posteriors occur for leaves with an ambiguous date such as `1986-XX-XX`, whose posterior is nearly flat inside the date range.

## Evidence

- Docstring diverges: the region should contain "the specified fraction of the probability mass" [packages/legacy/treetime/treetime/clock_tree.py#L1149-L1150](../../packages/legacy/treetime/treetime/clock_tree.py#L1149-L1150); the reported region holds all of it
- Objective not minimized: at the returned level the objective is `0.01`, while a level with objective `0` exists
- Adjacent correct: nodes whose posterior rises above height 10 on a side get a region with exactly 90% of the mass from the same code
- Reproduced: leaf `KP258696` of `rsv/a/20` (date `1986-XX-XX`, side heights at most `0.17` and `0.013`) in the golden-master capture `packages/treetime/src/timetree/__tests__/__fixtures__/gm_max_posterior_region_capture`. v0 reports `[1986.0014, 1986.9986]`, the whole support

## v0 impact

Nodes with a shallow posterior get the whole support as their 90% region, which overstates their date uncertainty. This includes `num_date_confidence` that augur writes from `get_max_posterior_region`.

## v1 status

The mass between the two level crossings is piecewise linear in `h`, so `fn highest_density_interval()` in [packages/treetime/src/timetree/max_posterior_region.rs](../../packages/treetime/src/timetree/max_posterior_region.rs) solves for the level exactly. For `KP258696`, v0's own CDF gives the v1 region 90% of the mass (`test_gm_max_posterior_region_shallow_posterior_holds_fraction_where_v0_returns_support`).
