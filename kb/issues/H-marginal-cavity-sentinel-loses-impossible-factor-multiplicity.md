# Marginal cavity sentinel loses impossible-factor multiplicity

Marginal inference represents an impossible aggregate with a scalar $-\infty$ sentinel. When multiple independent factors are impossible, removing one child can turn the aggregate into a finite neutral value even though another impossible factor remains.

`fn forward_log_lh_remove_child()` returns $0$ whenever both aggregate and removed-child scales are $-\infty$ [packages/treetime/src/partition/marginal_core.rs#L414-L426](../../packages/treetime/src/partition/marginal_core.rs#L414-L426). The scalar aggregate contains no count of other impossible factors.

## Counterexample

Let child factors $a$ and $b$ both contribute $-\infty$. Their combined log contribution is $-\infty$. A cavity operation that cancels only $a$ cannot infer whether $b$ remains from the scalar total, so it can incorrectly return zero.

## Options

- O1. Recompute the cavity aggregate while excluding the selected child. This directly preserves factor identity.
- O2. Store the finite aggregate plus the count of impossible factors. Removal decrements the count and remains impossible while the count is positive.

Recommendation: O2 if the representation can be made explicit and shared by dense/sparse code; O1 is the reference implementation for validation. The choice changes inference representation and requires approval before a ticket is executable.

## Recommendation

Adopt multiplicity-aware aggregation if it can be represented explicitly in both dense and sparse inference; use exclusion recomputation as the independent oracle. Obtain approval before creating an implementation ticket.

## Related issues

- [M-marginal-normalize-neg-infinity-masks-total.md](M-marginal-normalize-neg-infinity-masks-total.md)
- [M-sparse-marginal-zero-multiplicity-nan.md](M-sparse-marginal-zero-multiplicity-nan.md)
