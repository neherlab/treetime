# Marginal cavity sentinel loses impossible-factor multiplicity

> [!IMPORTANT]
> **Decision required.** The correction changes the inference representation and needs approval: recompute the cavity aggregate without the selected child (O1), or store the finite aggregate plus a count of impossible factors (O2). See "Options".

Marginal inference represents an impossible aggregate with a scalar $-\infty$ sentinel. When multiple independent factors are impossible, removing one child can turn the aggregate into a finite neutral value even though another impossible factor remains.

`fn forward_log_lh_remove_child()` returns `LogLh::ZERO` whenever both the aggregate and the removed child are `LogLh::IMPOSSIBLE` [packages/treetime/src/partition/marginal/shared/normalize.rs#L33-L39](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L33-L39). The scalar aggregate contains no count of other impossible factors. The dense forward pass ([packages/treetime/src/partition/marginal/shared/pass.rs#L404](../../packages/treetime/src/partition/marginal/shared/pass.rs#L404)) and the sparse forward pass ([packages/treetime/src/partition/marginal/sparse/forward.rs#L247](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L247)) both call it.

The defect needs at least two impossible child messages at one node. Positive branch lengths do not produce an impossible child message, so no reproduction on real input is known.

## Counterexample

Let child factors $a$ and $b$ both contribute $-\infty$. Their combined log contribution is $-\infty$. A cavity operation that cancels only $a$ cannot infer whether $b$ remains from the scalar total, so it can incorrectly return zero.

## Options

- O1. Recompute the cavity aggregate while excluding the selected child. This directly preserves factor identity.
- O2. Store the finite aggregate plus the count of impossible factors. Removal decrements the count and remains impossible while the count is positive.

Recommendation: O2 if the representation can be made explicit and shared by dense/sparse code; O1 is the reference implementation for validation. The choice changes inference representation and requires approval before implementation.

## Recommendation

Adopt multiplicity-aware aggregation if it can be represented explicitly in both dense and sparse inference; use exclusion recomputation as the independent oracle. Obtain approval before implementation starts.

## Related issues

- [M-marginal-normalize-neg-infinity-masks-total.md](M-marginal-normalize-neg-infinity-masks-total.md)
- [M-sparse-marginal-zero-multiplicity-nan.md](M-sparse-marginal-zero-multiplicity-nan.md)

## Proposed implementation plan

The locations above are current. The plan forks at the open decision.

1. Regression test first, in `packages/treetime/src/partition/marginal/shared/__tests__/test_marginal_shared.rs`: a node with two impossible child factors and one finite factor; removing one impossible child must stay `LogLh::IMPOSSIBLE` (oracle: hand calculation, the product of the remaining factors contains a zero). The existing case `test_marginal_shared_forward_log_lh_remove_child_cancels_matching_neg_infinity` asserts the single-factor cancel and stays valid only under O2 with count one
2. O1: in the dense forward pass (`packages/treetime/src/partition/marginal/shared/pass.rs#L404`) and the sparse forward pass (`packages/treetime/src/partition/marginal/sparse/forward.rs#L247`), compute the cavity log likelihood as the sum of the other children's message log likelihoods plus the parent's incoming message, and delete `fn forward_log_lh_remove_child()`
3. O2: replace the scalar `LogLh` of the node aggregate with a type that holds the finite sum and the count of impossible factors, built where the backward pass combines child messages; `fn forward_log_lh_remove_child()` decrements the count and returns `IMPOSSIBLE` while it stays positive. Use this type in both dense and sparse passes
4. Under O2, add a property test on random factor lists with random impossible entries: removal of any factor equals the recomputation without it (oracle: O1 recomputation)
5. Validate with `just lint-rs` and `just test-rs`; the marginal golden masters against v0 must not change, because no real input with positive branch lengths reaches the case
