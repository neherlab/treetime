# A resampled time distribution with two hard ends moves one of them by up to half a cell

The function-function convolution brings both operands to the finer of their two spacings (`fn on_common_spacing()`, [packages/treetime-distribution/src/distribution_ops/convolve.rs#L316-L326](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L316-L326)). `DistributionFunction::resample_dx()` ([packages/treetime-distribution/src/distribution_core/function.rs#L306-L314](../../packages/treetime-distribution/src/distribution_core/function.rs#L306-L313)) rounds the point count (`Grid::range_dx_point_count()`, [packages/treetime-grid/src/grid.rs#L155-L171](../../packages/treetime-grid/src/grid.rs#L155-L171)), keeps one end on a grid point and moves the other end by up to half a cell of the new spacing (`Grid::from_range_dx_anchored()`, [packages/treetime-grid/src/grid.rs#L100-L115](../../packages/treetime-grid/src/grid.rs#L100-L115)). `fn resample_anchor()` ([packages/treetime-distribution/src/distribution_core/function.rs#L425-L437](../../packages/treetime-distribution/src/distribution_core/function.rs#L425-L437)) keeps the single hard end of an operand with exactly one hard end, so branch distributions (hard at length 0, soft far end) and the usual subtree and parent distributions keep their hard end exactly ([kb/decisions/distribution-convolution-trapezoid-rule.md](../decisions/distribution-convolution-trapezoid-rule.md#operands-resampled-to-a-common-spacing-keep-their-single-hard-end-or-else-their-denser-end-on-a-grid-point)).

An operand with two hard ends keeps only its denser end (the left end when both ends agree), and its other hard end still moves by up to half a cell. The support of the convolution then ends up to half a cell before or after the true bound.

## Where it occurs

Release `timetree` with default flags (`--seed=42 --jobs=1`) on `flu/h3n2/500`, `ebola/362`, `zika/86`, `tb/149`, `rsv/a/2000`, `dengue/2000`, `lassa/L/500`, `mpox/clade-ii/2000` and `sc2/2844`, and the production cases `avian-ha`, `avian-pb2`, `oropouche-S` and `rubella-E1` (Nextstrain builds replayed with their own arguments, `--seed 1 --jobs 1`), counted with a temporary probe in `resample_dx()`. Only these resamplings moved a hard end:

| run | caller | operand | resamplings that moved a hard end | inward / outward | largest move (cells) | largest moved mass |
| --- | --- | --- | --- | --- | --- | --- |
| `avian-ha` | backward | subtree | 3773 | 949 / 2824 | 0.499 | 1.3e-3 |
| `avian-pb2` | backward | subtree | 800 | 290 / 510 | 0.499 | 2.1e-4 |
| `rubella-E1` | backward | subtree | 351 | 103 / 248 | 0.493 | 1.1e-1 |
| `mpox/clade-ii/2000` | forward | parent | 4 | 3 / 1 | 0.261 | 5e-123 |

The moved mass is the density at the moved end times the move, relative to the mass of the operand. The subtree operands are tips dated to a whole year (for example 2001.0014 to 2001.9986 in `rubella-E1`, on 3 or 4 grid points with a spacing of 0.33 years, resampled to 0.115 years), whose time distribution is the flat date range with two hard ends. No branch distribution with two hard ends (an undecayed likelihood that ends `Hard` at `MAX_BRANCH_LENGTH`) was resampled in these runs.

In a backward message the output window clips the support at the bound of the unresampled operands (`WindowedConvolution::new()`, [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L259-L268](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L259-L268)), so an outward move of the right end has no effect on the end of the message, and an inward move makes the message end before the latest date that the parent can take.

The effect of these moves on the node dates is not measured alone. Keeping the single hard end changed the dates of every run above except `tb/149` (largest leaf and matched internal node shifts of 0.7 days on `zika/86` to about 2000 days on `rubella-E1`, with tree changes on 9 of the 13 runs), which mixes this remainder with the moves that the single-hard-end rule removed.

> [!IMPORTANT]
> **Decision required.** The remaining moves come from operands with two hard ends. Options:
>
> - Keep both ends exact: choose the common spacing as range/⌈range/dx⌉ of the operand with two hard ends and resample the finer operand to it, a change of `on_common_spacing()` in `convolve.rs`. Both supports then stay exact, at a spacing finer than the finer of the two spacings by less than one part in the cell count of that operand
> - Accept the half-cell error at the second hard end, and document it in [kb/decisions/distribution-convolution-trapezoid-rule.md](../decisions/distribution-convolution-trapezoid-rule.md)
>
> Evidence: the table above. The moved mass reaches 11% of a year-only tip date range in `rubella-E1`, where the range covers only 3 or 4 grid points.

## Related issues

- [H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md](H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md): the accuracy of backward messages along zero-length branches
