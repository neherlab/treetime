# Branch distribution grid uses uniform spacing

The branch distribution grid in `fn create_simple_grid()` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L69-L80](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L69-L80) uses `Array1::linspace` (uniform spacing) over `[min_bl, min(max(center * 5, one_mutation * 10), MAX_BRANCH_LENGTH)]`, where `min_bl` is 0 when the likelihood has a finite value at zero length and `one_mutation * 0.01` otherwise. The grid is now branch-length-informed (extent scales with the ML branch length, capped at 5.0 subs/site) and uses 300 points, a substantial improvement over the prior 1000-point grid that spread resolution across a 200-year time window.

## Current status

Mass sizing now governs the final extent: `create_simple_grid()` builds the uniform pilot grid described above, then `fn rewindow_to_mass()` resamples it onto a mass-bounded domain [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L63](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L63). The `min(center * 5, ...)` extent above is therefore the pilot grid the rewindow measures, not the stored grid. The residual concern is unchanged: the pilot grid is uniform `linspace`, so a too-narrow or under-resolved pilot can clip or blur the peak before the rewindow measures mass. This clipping (mechanism (a) in [H-timetree-mass-sizing-node-times-break-downstream-invariants.md](H-timetree-mass-sizing-node-times-break-downstream-invariants.md)) was investigated there and found not to be the cause of the observed `+inf`, which is a topology inversion; the clipping remains a latent resolution and accuracy risk, not a confirmed failure.

The resolution improvement has not been empirically validated: the golden master tests (`test_gm_runner_marginal_dense`, `test_gm_runner_marginal_sparse`, `test_gm_runner_poisson`) remain `#[ignore]`d at target tolerance `epsilon = 1e-6`. Un-ignoring them and verifying that they pass at that tolerance is the acceptance criterion for closing this issue.

> [!WARNING]
> **Needs review.** The `#[ignore]` reasons of `test_gm_runner_marginal_dense` [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L45](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L45) and `test_gm_runner_marginal_sparse_dense_consistency` [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L45](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L45) state a gap of 0.92 years. A stack of 2026-10-18 measured the dense gap to v0 on `flu_h3n2_20` as 0.057376 years at `NODE_0000017`, with both indel policies, after it added an exact mass model and an exact convolution of ranges to `rust` of that day and before it changed the convolution grid; neither change is on `rust`. The gap on `rust` of 2026-10-19 has not been measured; run the dense test with `--run-ignored only` before relying on either number.

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

Both options need either non-uniform grid support in `DistributionFunction` (the likelihood is now built through `GridFn::from_range_values()`, which takes a range and therefore a uniform grid) or a resample step onto a uniform grid before the distribution is built.

> [!IMPORTANT]
> **Decision required.** The choice between the two grid designs is open, and it applies only if the golden master tests still fail once the grid-width discrepancy is the only remaining cause. Option one, v0 segments, keeps parity with v0's `BranchLenInterpolator` grid and makes the grid easier to compare against v0. Option two, two tiers, is simpler and keeps most of the current uniform machinery, but it diverges from v0 and needs an entry in `kb/decisions/`. Evidence: the tests are currently ignored for discrepancies of 0.27 years ([packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L35](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L35)) and 0.92 years at a root-adjacent node ([packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L45](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L45)).
