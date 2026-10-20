# Branch distribution grid uses uniform spacing

The branch distribution grid in `fn create_simple_grid()` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L69-L80](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L69-L80) uses `Array1::linspace` (uniform spacing) over `[min_bl, min(max(center * 5, one_mutation * 10), MAX_BRANCH_LENGTH)]`, where `min_bl` is 0 when the likelihood has a finite value at zero length and `one_mutation * 0.01` otherwise. The extent scales with the ML branch length, is capped at 5.0 subs/site, and holds 300 points.

## Current status

Mass sizing governs the final extent: `create_simple_grid()` builds the uniform pilot grid described above, then `fn rewindow_to_mass()` resamples it onto a mass-bounded domain [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L63](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L63). The `min(center * 5, ...)` extent above is therefore the pilot grid the rewindow measures, not the stored grid. The remaining concern: the pilot grid is uniform `linspace`, so a too-narrow or under-resolved pilot can clip or blur the peak before the rewindow measures mass. This clipping does not cause the `+inf` values that follow from inverted node times; those belong to [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md). The clipping remains a latent resolution and accuracy risk, not a confirmed failure.

## Measured gap to v0

`test_gm_runner_marginal_dense` [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L30-L32](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L30-L32) runs at `epsilon = 1e-1` and cites this issue: v1 dense node times differ from v0 by up to 0.070 years on `flu_h3n2_20` (`NODE_0000017`, both indel policies). The gaps of the other datasets are in [N-timetree-gm-runner-dataset-cases-disabled.md](N-timetree-gm-runner-dataset-cases-disabled.md).

> [!IMPORTANT]
> **Investigation required.** No measurement shows that the uniform pilot grid causes these gaps. Isolate the cause before the grid decision below, for example by running the dense test with a finer or a v0-shaped pilot grid and comparing the gap.

The acceptance criterion for closing this issue: `test_gm_runner_marginal_dense` passes at `epsilon = 1e-6`.

## Remaining concern

For very short branches where `one_mutation * 10` governs the extent, the uniform grid may still under-resolve the peak relative to v0's non-uniform grid. A non-uniform grid concentrating points near the peak while maintaining broad tail coverage (as v0 does) would further improve resolution for these cases.

## v0 reference

v0 uses a 5-segment non-uniform grid (`branch_len_interpolator.py:50-62`):

- Log-spaced near zero (5 points)
- Linear near zero (8 points)
- Quadratic from 0 to peak (n/3 points, dense near peak)
- Quadratic from peak to 3\*sigma (n/3 points)
- Quadratic from 3\*sigma to MAX_BRANCH_LENGTH (n/3 points, sparse tail)

This concentrates ~40 of 125 points near the peak while covering up to MAX_BRANCH_LENGTH = 4.0 subs/site.

## Fix options

If the golden master tests pass at `epsilon = 1e-6` with the uniform pilot grid, the uniform grid is sufficient and the issue closes without a non-uniform grid. Otherwise the pilot grid needs more resolution at the peak:

- **v0 segments**: a non-uniform grid matching v0's multi-segment design above
- **Two tiers**: a fine uniform grid over the peak region, concatenated with a coarser grid for the tail

Both options need either non-uniform grid support in `DistributionFunction` (the likelihood is built through `GridFn::from_range_values()`, which takes a range and therefore a uniform grid) or a resample step onto a uniform grid before the distribution is built.

> [!IMPORTANT]
> **Decision required.** The choice between the two grid designs is open, and it applies only if the golden master tests still fail once the grid-width discrepancy is the only remaining cause. Option one, v0 segments, keeps parity with v0's `BranchLenInterpolator` grid and makes the grid easier to compare against v0. Option two, two tiers, is simpler and keeps most of the current uniform machinery, but it diverges from v0 and needs an entry in `kb/decisions/`. Evidence: the tests are ignored for gaps to v0 of 0.070 years (dense) and up to 0.349 years (Poisson) on `flu_h3n2_20`, see Measured gap to v0. No measurement yet attributes these gaps to the grid.
