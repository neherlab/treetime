# Marginal forward pass divides in probability space with a zero-divisor floor

> [!IMPORTANT]
> **Investigation required.** No wrong result is shown. Reproduce the reachable zero patterns of numerator and divisor together, derive the intended cavity-message semantics, and compare the dense forward profiles with a log-space recomputation, as listed under "Investigation required".

The forward pass forms the message to a child by dividing the node posterior by the message from that child. The dense, discrete, and sparse passes do this in plain probability space, through one shared helper: `fn divide_out_value()` returns `value / divisor.max(f64::MIN_POSITIVE)` [packages/treetime/src/partition/marginal/shared/normalize.rs#L51-L53](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L51-L53). The dense `fn divide_out()` [packages/treetime/src/partition/marginal/shared/normalize.rs#L45-L49](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L45-L49) and the sparse `fn divide_out_site()` [packages/treetime/src/partition/marginal/sparse/message.rs#L184-L188](../../packages/treetime/src/partition/marginal/sparse/message.rs#L184-L188) call it. The dense forward pass divides and then normalizes [packages/treetime/src/partition/marginal/shared/pass.rs#L320-L321](../../packages/treetime/src/partition/marginal/shared/pass.rs#L320-L321).

v0's preorder divides in log space and multiplies back in probability space [packages/legacy/treetime/treetime/treeanc.py#L880-L917](../../packages/legacy/treetime/treetime/treeanc.py#L880-L917). Both compute the same quantity, but the floating-point paths differ.

## Floor of structural zeros

When the message from a child has an exact zero (a state is impossible given the subtree evidence), plain division gives infinity, so the floor replaces the divisor with `f64::MIN_POSITIVE` (about 2.2e-308). The claimed support change is not established: an exact-zero numerator divided by `f64::MIN_POSITIVE` remains zero, while a positive numerator over an exact-zero divisor indicates a different, inconsistent-message case. The numerator and divisor states must be reproduced together before a failure mechanism can be assigned.

## Plain probability space

Division by near-zero probabilities can amplify rounding errors at positions where a child gives near-zero likelihood to some states. Two facts limit the concern:

- The backward pass already combines child messages in log space (`normalize_from_log` in the dense pass [packages/treetime/src/partition/marginal/shared/pass.rs#L176-L180](../../packages/treetime/src/partition/marginal/shared/pass.rs#L176-L180), `softmax_with_log_norm_inplace()` in the sparse pass), which removes the polytomy underflow mode that a public v0 report describes [[issue](https://github.com/neherlab/treetime/issues/39)]; the maintainer named log-space message accumulation as the remedy [[comment](https://github.com/neherlab/treetime/issues/39#issuecomment-383567696)]
- The parent profile already contains the child message as a factor, so a near-zero divisor pairs with a near-zero numerator; only exact zeros reach the floor

The golden-master tests compare v1 dense marginal reconstruction with v0 at tolerances of 1e-6 to 1e-7, which absorbs the differences of the two paths.

## Impact

No output change is shown. The mass that the floor can add is below 1e-300. The concern is principled: the operation can change the support of a probability vector, and plain division can lose precision where log-space subtraction would not.

## Investigation required

- Reproduce exact-zero and subnormal divisor cases with their corresponding posterior numerators
- Derive the cavity-message limit for each reachable zero pattern
- Compare the dense forward profiles with a log-space recomputation on a large polytomy and on extreme branch lengths
- Compare the dense, sparse, and v0 log-space paths against an independent analytical case

## Readiness

No correction is ready until the reachable zero patterns and the intended cavity-message semantics are established.

## Related issues

- [M-marginal-cavity-sentinel-loses-impossible-factor-multiplicity.md](M-marginal-cavity-sentinel-loses-impossible-factor-multiplicity.md)
- [M-sparse-marginal-zero-multiplicity-nan.md](M-sparse-marginal-zero-multiplicity-nan.md)
- [N-distribution-plain-division-fixed-floor.md](N-distribution-plain-division-fixed-floor.md)
