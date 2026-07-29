# Marginal forward pass zero-divisor floor converts structural zeros

The forward pass msg_to_child computation divides the node posterior by the child's msg_from_child. When msg_from_child contains exact zeros (a state is structurally impossible given the subtree evidence), the division would produce infinity. The code floors the divisor to `f64::MIN_POSITIVE` (~2.2e-308) to avoid this.

The claimed support change is not established: an exact-zero numerator divided by `f64::MIN_POSITIVE` remains zero, while a positive numerator over an exact-zero denominator indicates a different inconsistent-message case. The numerator and divisor states must be reproduced together before assigning a failure mechanism.

## Affected code

- Sparse: [packages/treetime/src/partition/marginal_passes.rs#L295](../../packages/treetime/src/partition/marginal_passes.rs#L295)
- Dense + discrete (shared): [packages/treetime/src/partition/marginal_core.rs#L163](../../packages/treetime/src/partition/marginal_core.rs#L163)

## Impact

The practical impact is negligible because `f64::MIN_POSITIVE` is ~2.2e-308, making the resulting probability mass undetectable by downstream computations. The concern is principled: the operation changes the topology of the probability simplex (zero support becomes nonzero support).

## Investigation required

- Reproduce exact-zero and subnormal divisor cases with their corresponding posterior numerators.
- Derive the cavity-message limit for each reachable zero pattern.
- Compare dense, sparse, and v0 log-space paths against an independent analytical case.

## Ticket readiness

No correction is ready until the reachable zero patterns and intended cavity-message semantics are established.
