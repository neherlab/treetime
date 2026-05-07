# Mugration optimize_gtr_rate restores mu without restoring profiles

`optimize_gtr_rate` at [packages/treetime/src/commands/mugration/gtr_refinement.rs#L198-L218](../../packages/treetime/src/commands/mugration/gtr_refinement.rs#L198-L218) evaluates three cost points (`cost_lo`, `cost_mid`, `cost_hi`) by running `run_backward`, which mutates the partition (sets `mu`, clears `log_lh`, runs backward pass to update all node/edge profiles). When no interior minimum exists, the function restores `old_mu` at line 216 but leaves the backward-pass profiles in their "evaluated at `hi`" state (where `hi = 100 * sqrt(old_mu)`).

The next call inside `refine_gtr_iterative` is `count_transitions_discrete`, which reads edge `msg_to_parent`/`msg_to_child`/`msg_from_child`. These profiles were produced under a different `mu` than the restored `old_mu`, making the transition counts internally inconsistent: `expQt(branch_length)` uses the restored `mu`, but the message profiles it multiplies were computed under the `hi` value.

## v0 comparison

v0 reproduces the same inconsistency (scipy Brent evaluates the cost function at multiple points, leaving profiles in the state of the last evaluation). This is not a v0/v1 divergence but rather a shared correctness issue.

## Affected code

- `mu` restore without profile restore: [packages/treetime/src/commands/mugration/gtr_refinement.rs#L213-L218](../../packages/treetime/src/commands/mugration/gtr_refinement.rs#L213-L218)
- Profile consumers: `count_transitions_discrete` at [gtr_refinement.rs#L121-L125](../../packages/treetime/src/commands/mugration/gtr_refinement.rs#L121-L125)

## Fix

After restoring `old_mu`, re-run `run_backward` to recompute profiles consistent with the restored `mu`. The cost of one extra backward pass per no-bracket iteration is small relative to the correctness guarantee.
