# Mugration optimize_gtr_rate restores mu without restoring profiles

> [!IMPORTANT]
> **Decision required.** The current behavior is v0 parity, and the fix below diverges from v0. v0's `optimize_gtr_rate()` restores `old_mu` after a failed bracket but keeps the profiles of the last cost evaluation ([packages/legacy/treetime/treetime/treeanc.py#L1679-L1708](../../packages/legacy/treetime/treetime/treeanc.py#L1679-L1708)). The mugration golden masters match v0 within 1e-7 with the current behavior ([kb/algo/mugration.md](../algo/mugration.md#golden-master-tolerance)); a change here moves the results away from v0 where a bracket fails. The options are:
>
> - Keep v0 parity and record the inconsistency as v0 behavior
> - Recompute the backward pass at the restored `mu`, as described under Fix, and record the divergence from v0

`fn optimize_gtr_rate()` [packages/treetime/src/gtr/refinement.rs#L156-L264](../../packages/treetime/src/gtr/refinement.rs#L156-L264) evaluates three cost points (`cost_lo`, `cost_mid`, `cost_hi`). Each evaluation runs the marginal backward pass under a candidate `mu` and returns the candidate node states and backward messages. Its parameter `range_edge` selects what happens when no interior minimum exists. With `RateRangeEdge::KeepStart`, which mugration uses, `fn keep_start_rate()` returns the candidate evaluated at `hi = 100 * sqrt(old_mu)` with its `mu` reset to `old_mu` ([packages/treetime/src/gtr/refinement.rs#L272-L292](../../packages/treetime/src/gtr/refinement.rs#L272-L292)). Only when the `hi` evaluation failed does it run a new backward pass at `old_mu`.

The next iteration of `fn refine_gtr()` calls `count_transitions` with these backward messages ([packages/treetime/src/gtr/refinement.rs#L45-L53](../../packages/treetime/src/gtr/refinement.rs#L45-L53)), and the next rate step starts from these node states ([packages/treetime/src/gtr/refinement.rs#L63](../../packages/treetime/src/gtr/refinement.rs#L63)). The messages were produced under the `hi` rate, while `expQt(branch_length)` uses the restored `mu`, so the transition counts are internally inconsistent.

A related mix of models affects every iteration, also those with a bracket: `fn refine_gtr()` with `RefinementRate::Optimize` runs one full marginal update under the first fitted model and passes its `forward` messages and root state to `count_transitions` in every later iteration ([packages/treetime/src/gtr/refinement.rs#L44-L75](../../packages/treetime/src/gtr/refinement.rs#L44-L75)). The loop replaces the backward messages after each rate step but does not recompute the forward messages, so from the second fit on they belong to the first fitted model and rate. This is the v0 order ([kb/algo/mugration.md](../algo/mugration.md#gtr-model-construction)); a full update per iteration is proposed in [kb/proposals/mugration-full-reconstruction-per-iteration.md](../proposals/mugration-full-reconstruction-per-iteration.md).

The relative-rate step of optimize and timetree with several alignment files uses `RateRangeEdge::LowerCostEnd` instead: it returns the end of the range with the lower cost together with the node states and messages computed at that rate, and the caller runs a new marginal update with the returned rates ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)). That step has no inconsistency. Mugration keeps `RateRangeEdge::KeepStart` until the decision above is made.

## v0 comparison

v0 reproduces the same inconsistency: scipy Brent evaluates the cost function at the bracket points, and the bracket failure leaves the profiles in the state of the last evaluation. The inconsistency is shared by v0 and v1, not a divergence between them.

## Fix

After restoring `old_mu`, run the backward pass again to recompute profiles consistent with the restored `mu`. The cost of one extra backward pass per no-bracket iteration is small relative to the correctness guarantee.
