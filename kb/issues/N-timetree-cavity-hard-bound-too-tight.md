# Cavity division cannot restore a hard bound the product already removed

The forward-pass cavity is computed as `parent_posterior / msg_to_parent`: `fn refine_distribution_from_parent()` calls `distribution_division()` with the parent posterior and the message of the child ([packages/treetime/src/timetree/inference/forward_pass.rs#L193-L196](../../packages/treetime/src/timetree/inference/forward_pass.rs#L193-L196)). The posterior is the product of all the parent's factors, and the message is one of them, so `divide(f, g)` recovers the product of the remaining factors. On the hard sides this recovery is exact only when the divided-out child does **not** hold the tightest hard bound.

## Symptom and scope

When the divided-out child `g` held the tightest hard bound on a side, the parent product `f` stored only that child's bound; it never stored the next factor's (looser) bound. Dividing `f` by `g` is then a `0/0` cancellation at that edge, and the quotient can only keep `f`'s hard bound there -- the former child's bound -- rather than the looser bound the remaining factors actually imply. The cavity is therefore slightly too tight on that side: it excludes a sliver of support the rest of the tree still permits.

This is a bounded inexactness, not a collapse. It affects only the hard side whose tightest bound belonged to the removed child, and only by the gap between that child's bound and the next-tightest factor's bound. The soft-tail fix (`fn divide_function_by_function()`, see [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md)) resolves the empty-distribution defect regardless; this residual tightness is a separate, orthogonal concern.

## Root cause

Division operates on the aggregate posterior, which has already discarded the per-factor hard bounds during multiplication. A pure division cannot reconstruct information the product did not retain.

## Fix approach

Rebuild the cavity by multiplying all factors _except_ the selected child, instead of dividing the aggregate. That restores the exact hard support and avoids the `0/0` cancellation entirely.

The backward pass keeps the message of every child: `NodeBackward.to_parent` holds the message of each node to its parent ([packages/treetime/src/timetree/inference/result.rs#L14-L18](../../packages/treetime/src/timetree/inference/result.rs#L14-L18)). The "all but one" product also needs the date constraint of the parent and the message that the parent received from its own parent. `NodePosterior` does not store that message; the forward pass multiplies it into the parent posterior and drops it ([packages/treetime/src/timetree/inference/forward_pass.rs#L198-L209](../../packages/treetime/src/timetree/inference/forward_pass.rs#L198-L209)).

This is a forward-pass structural change, out of scope for the division operator itself.

> [!WARNING]
> **Needs review.** Since the commit `refactor(timetree): run the time passes on tree scans and free posterior grids early`, the forward tree scan hands the backward data of each node only to the forward step of that node ([packages/treetime/src/timetree/inference/forward_pass.rs#L38-L42](../../packages/treetime/src/timetree/inference/forward_pass.rs#L38-L42), [packages/treetime/src/timetree/inference/forward_pass.rs#L98](../../packages/treetime/src/timetree/inference/forward_pass.rs#L98)). The forward step of a child therefore does not see the messages of its siblings, which the "all but one" product needs. Check how the step can receive them before acting on the paragraph above.

> [!IMPORTANT]
> **Decision required.** The exact cavity costs grid multiplications and memory that the division does not, and the division matches v0. Options:
>
> - **Keep the division**: v0 parity, no extra cost. The cavity stays slightly too tight on the hard side whose tightest bound belonged to the removed child
> - **Rebuild the cavity from factors**: keep the message from the parent of each internal node, and multiply it with the date constraint and the messages of the other children. For a node with `k` children, each cavity costs `k` multiplications instead of one division; reusing prefix and suffix products of the child messages reduces this to about three multiplications per child
>
> Evidence: the division at [packages/treetime/src/timetree/inference/forward_pass.rs#L193-L196](../../packages/treetime/src/timetree/inference/forward_pass.rs#L193-L196); the stored child messages above; the v0 contract and the prototype test below

## Reference

- v0 defines division only for a numerator that contains the denominator (`packages/legacy/treetime/treetime/distribution.py` `Distribution.divide()`), so the cavity contract (dividend contains divisor) is the documented one; v0 has the same limitation.
- The Python prototype `test_scripts/density_algebra.py` `test_dividing_out_the_boundary_factor_cannot_restore_support` demonstrates the limitation directly.
