# Floating-point assertion tolerance is hidden by a variable

`assert_sparse_profile_normalized` accepts `max_ulps` and forwards it to floating-point assertions. Current callers pass a literal, but the assertion site does not expose or constrain the oracle tolerance, violating the repository’s test policy.

The helper parameter is forwarded into `pretty_assert_ulps_eq!` in [packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs#L35-L59](../../packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs#L35-L59).

## Potential solutions

- O1. Put explicit literal tolerances at every assertion.
- O2. Define one fixed domain assertion whose tolerance and oracle derivation are part of the helper contract.

## Recommendation

Use explicit literal `max_ulps` values at every assertion site. Remove the variable tolerance parameter; do not introduce a helper that hides the literal or widen the value to make failures pass.

## Related issues

- [N-test-quality-deficiencies.md](N-test-quality-deficiencies.md)
- [N-test-coverage-gaps.md](N-test-coverage-gaps.md)
