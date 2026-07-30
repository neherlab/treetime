# Backward pass child-order permutation tests missing

The backward pass multiplies child parent-time messages sequentially. The result should be independent of processing order (multiplication is commutative). No test currently verifies this invariant.

A property test should:

- Generate N child messages with random supports and Constant left / Zero right tails
- Multiply them in multiple permutations with normalize between steps
- Assert the final distributions have the same `likely_time()` across all orderings

This would catch ordering-dependent failures from tail metadata loss, normalize interactions, or asymmetric grid resolution.

## Related

- [kb/issues/M-distribution-multiply-normalize-tail-composition.md](M-distribution-multiply-normalize-tail-composition.md): the underlying tail composition issue that makes ordering matter
