# get_max_posterior_region mixes modes on a non-monotone side of the peak

## v0 location

`ClockTree.get_max_posterior_region()` (`#get_max_posterior_region`) [packages/legacy/treetime/treetime/clock_tree.py#L1183-L1211](../../packages/legacy/treetime/treetime/clock_tree.py#L1183-L1211)

## Erratum

For a node whose marginal time posterior has its peak inside the grid, v0 finds the 90% region by mapping a level `h` above the peak to a position on each side of the peak. It builds this map with `scipy.interpolate.interp1d(heights, positions)` from all points of that side. `interp1d` sorts the points by height. When a side is not monotone, for example a second mode or a shoulder, neighbors in height order come from different parts of the side. The interpolated boundary can then lie beyond a dip in the density, and the region contains points whose density is below the level of its own boundaries.

Example: a side walked outward from the peak with heights `0, 2, 1, 3` at positions `0, 1, 2, 3`. Level 1.5 maps to position 1.5. The first crossing of level 1.5 outward from the peak is at position 0.75. The region up to 1.5 contains position 1, at height 2, which is above the level.

## Evidence

- Comment diverges: the code states its own assumption, "this assumes there is only one peak --- might fail in odd cases" [packages/legacy/treetime/treetime/clock_tree.py#L1185](../../packages/legacy/treetime/treetime/clock_tree.py#L1185)
- Docstring diverges: the docstring describes "the interval around the highest posterior probability region" [packages/legacy/treetime/treetime/clock_tree.py#L1149](../../packages/legacy/treetime/treetime/clock_tree.py#L1149). A region that contains points below its boundary level is not a highest-density region
- Reproduced: the example above, run with the scipy of the v0 Python environment (`./dev/docker/python`)

## v0 impact

Only posteriors with a non-monotone side are affected. On every single-peaked posterior the sorted-height map equals the first crossing outward from the peak. Among the v0 posteriors of `ebola/20` and `rsv/a/20`, one internal node has a non-monotone side, at height 13.8 above the peak, far above the 90% level. Its region is the same in v0 and v1.

## v1 status

`fn max_posterior_region()` in [packages/treetime/src/timetree/max_posterior_region.rs](../../packages/treetime/src/timetree/max_posterior_region.rs) replaces the heights of each side by their running maximum outward from the peak. A level then maps to its first crossing outward from the peak, so the region is the level set connected to the peak. A level above a dip adds the mode beyond it as a whole. On monotone sides the result equals v0, except for the mass of tail laws and hard cells, which v1 counts ([kb/decisions/timetree-posterior-statistics-count-tail-mass.md](../decisions/timetree-posterior-statistics-count-tail-mass.md)).
