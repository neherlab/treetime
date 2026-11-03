# A backward time message can miss the date of its child by up to half a cell

The backward pass sends the time message of a node to its parent by convolving the time distribution of the node given its subtree with the negated branch-time distribution (`fn send_backward_message()`, [packages/treetime/src/timetree/inference/backward_pass.rs#L165-L193](../../packages/treetime/src/timetree/inference/backward_pass.rs#L165-L193)). The negated branch distribution has a hard right end at branch length 0, so the message should end exactly where the subtree distribution ends: a parent can be as late as its child, but not later.

The function-function convolution first brings both operands to the finer of their two spacings (`fn on_common_spacing()`, [packages/treetime-distribution/src/distribution_ops/convolve.rs#L316-L326](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L316-L326)). When the branch distribution has the coarser spacing, `DistributionFunction::resample_dx()` ([packages/treetime-distribution/src/distribution_core/function.rs#L306-L323](../../packages/treetime-distribution/src/distribution_core/function.rs#L306-L323)) resamples it on a grid anchored at its denser end, and at the left end when both ends have the same value. The point count is rounded (`Grid::range_dx_point_count()`, [packages/treetime-grid/src/grid.rs#L155-L171](../../packages/treetime-grid/src/grid.rs#L155-L171)), so the end opposite the anchor moves by up to half a cell (`Grid::from_range_dx_anchored()`, [packages/treetime-grid/src/grid.rs#L100-L115](../../packages/treetime-grid/src/grid.rs#L100-L115)). When the far end of the negated branch distribution (the longest branch length) has the higher density, or both ends agree, the hard end at 0 is the end that moves, and the message ends up to half a cell of the common spacing before or after the end of the subtree distribution.

Commit `8fa5692db` introduced the anchor at the denser end for the case where the zero end is the denser one: a short branch has its highest density at length 0, and the resampled grid must keep that end on a grid point. The rule does not cover a hard end that is not the denser end.

## Example

In a replayed Nextstrain avian influenza HA build, the subtree distribution of an internal node ends at 2024.6926229508197, its branch distribution ends at 0, and the backward message ends at 2024.69196: 0.00066 years (about 6 hours) before the latest date the parent can take. The issue exists on `rust` before and after the change that keeps both grid ends as given coordinates.

> [!IMPORTANT]
> **Investigation required.** The effect on node dates and on the log likelihood of the node times is not measured. Collect: how often the hard end moves in the backward passes of the example datasets and of the production case above, in which direction (a message that ends after the date of the child admits a parent later than its child), the probability mass the message gains or loses in the moved part, and the change of the node dates when the hard end stays on a grid point.

## Fix direction

Anchor a resampled operand at its hard end when it has one, and at its denser end otherwise. For the negated branch distribution of a short branch, both rules choose the end at 0, so the error that `8fa5692db` fixed stays fixed. Add a convolution test with a branch distribution whose far end is denser and whose spacing is coarser than that of the subtree distribution, and assert that the message ends exactly at the end of the subtree distribution.

## Related issues

- [H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md](H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md): the accuracy of backward messages along zero-length branches
