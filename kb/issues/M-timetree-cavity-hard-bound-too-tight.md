# Cavity division cannot restore a hard bound the product already removed

The forward-pass cavity is computed as `parent_posterior / msg_to_parent` (`fn refine_node()` in `packages/treetime/src/timetree/inference/forward_pass.rs`). The posterior is the product of all the parent's factors, and the message is one of them, so `divide(f, g)` recovers the product of the remaining factors. On the hard sides this recovery is exact only when the divided-out child does **not** hold the tightest hard bound.

## Symptom and scope

When the divided-out child `g` held the tightest hard bound on a side, the parent product `f` stored only that child's bound; it never stored the next factor's (looser) bound. Dividing `f` by `g` is then a `0/0` cancellation at that edge, and the quotient can only keep `f`'s hard bound there -- the former child's bound -- rather than the looser bound the remaining factors actually imply. The cavity is therefore slightly too tight on that side: it excludes a sliver of support the rest of the tree still permits.

This is a bounded inexactness, not a collapse. It affects only the hard side whose tightest bound belonged to the removed child, and only by the gap between that child's bound and the next-tightest factor's bound. The soft-tail fix (`fn divide_function_by_function()`, see [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md)) resolves the empty-distribution defect regardless; this residual tightness is a separate, orthogonal concern.

## Root cause

Division operates on the aggregate posterior, which has already discarded the per-factor hard bounds during multiplication. A pure division cannot reconstruct information the product did not retain.

## Fix approach

Rebuild the cavity by multiplying all factors _except_ the selected child, instead of dividing the aggregate. This requires the forward pass to retain each child's contribution (the per-child factors) rather than only the combined posterior, so the "all but one" product can be reformed on demand. That restores the exact hard support and avoids the `0/0` cancellation entirely.

This is a forward-pass structural change (factor retention), out of scope for the division operator itself. Confirmed as the correct exact approach in prior analysis; deferred pending a decision on the memory cost of retaining per-child factors versus recomputing the product.

## Reference

- v0 defines division only for a numerator that contains the denominator (`packages/legacy/treetime/treetime/distribution.py` `Distribution.divide()`), so the cavity contract (dividend contains divisor) is the documented one; v0 has the same limitation.
- The Python prototype `test_scripts/density_algebra.py` `test_dividing_out_the_boundary_factor_cannot_restore_support` demonstrates the limitation directly.
