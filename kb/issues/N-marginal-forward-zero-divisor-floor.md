# Marginal forward pass zero-divisor floor converts structural zeros

The forward pass msg_to_child computation divides the node posterior by the child's msg_from_child. When msg_from_child contains exact zeros (a state is structurally impossible given the subtree evidence), the division would produce infinity. The code floors the divisor to `f64::MIN_POSITIVE` (~2.2e-308) to avoid this.

The floor converts structurally impossible states (probability exactly zero) into states with tiny but nonzero support in the msg_to_child. This changes the information sent to the child: instead of "this state is impossible given the outgroup", the message becomes "this state has negligible but nonzero probability".

## Affected code

- Sparse: [packages/treetime/src/representation/partition/marginal_passes.rs#L295](../../packages/treetime/src/representation/partition/marginal_passes.rs#L295)
- Dense: [packages/treetime/src/representation/partition/marginal_dense.rs#L324](../../packages/treetime/src/representation/partition/marginal_dense.rs#L324)
- Discrete: [packages/treetime/src/representation/partition/discrete.rs#L206](../../packages/treetime/src/representation/partition/discrete.rs#L206)

## Impact

The practical impact is negligible because `f64::MIN_POSITIVE` is ~2.2e-308, making the resulting probability mass undetectable by downstream computations. The concern is principled: the operation changes the topology of the probability simplex (zero support becomes nonzero support).

## Alternative approaches

- Log-space division (subtraction) with explicit handling of -inf - -inf
- Propagate structural zeros through the forward pass by skipping positions where msg_from_child is zero
- v0 uses log-space arithmetic in the forward pass, avoiding the division entirely
