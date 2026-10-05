# Mugration optimize_gtr_rate restores mu without restoring profiles

> [!IMPORTANT]
> **Decision required.** The current behavior is v0 parity, and the fix below diverges from v0. v0's `optimize_gtr_rate()` restores `old_mu` after a failed bracket but keeps the profiles of the last cost evaluation ([packages/legacy/treetime/treetime/treeanc.py#L1679-L1708](../../packages/legacy/treetime/treetime/treeanc.py#L1679-L1708)). The rate-optimizer evaluation point is also a candidate cause of the open v0 parity gap in [M-mugration-iterative-gtr.md](M-mugration-iterative-gtr.md), so a change here interacts with that issue. The options are:
>
> - Keep v0 parity and record the inconsistency as v0 behavior
> - Recompute the backward pass at the restored `mu`, as described under Fix, and record the divergence from v0

`fn optimize_gtr_rate()` [packages/treetime/src/gtr/refinement.rs#L129-L232](../../packages/treetime/src/gtr/refinement.rs#L129-L232) evaluates three cost points (`cost_lo`, `cost_mid`, `cost_hi`). Each evaluation runs the marginal backward pass under a candidate `mu` and returns the candidate node states and backward messages. When no interior minimum exists, the function returns the candidate evaluated at `hi = 100 * sqrt(old_mu)` with its `mu` reset to `old_mu` ([packages/treetime/src/gtr/refinement.rs#L219-L231](../../packages/treetime/src/gtr/refinement.rs#L219-L231)). Only when the `hi` evaluation failed does it run a new backward pass at `old_mu`.

The next iteration of `fn refine_gtr_model_and_rate()` calls `count_transitions` with these node states and backward messages ([packages/treetime/src/gtr/refinement.rs#L79-L83](../../packages/treetime/src/gtr/refinement.rs#L79-L83)). The messages were produced under the `hi` rate, while `expQt(branch_length)` uses the restored `mu`, so the transition counts are internally inconsistent.

## v0 comparison

v0 reproduces the same inconsistency: scipy Brent evaluates the cost function at the bracket points, and the bracket failure leaves the profiles in the state of the last evaluation. The inconsistency is shared by v0 and v1, not a divergence between them.

## Fix

After restoring `old_mu`, run the backward pass again to recompute profiles consistent with the restored `mu`. The cost of one extra backward pass per no-bracket iteration is small relative to the correctness guarantee.
