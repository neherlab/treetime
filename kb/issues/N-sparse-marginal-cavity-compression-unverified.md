# Sparse cavity-profile compression lacks an error contract

## Summary

`fn compute_msg_to_child()` stores normalized cavity profiles explicitly even when their MAP state matches the fixed-row fallback and their probability is near $1$. Compressing those profiles approximately could improve representation density, but no analytical error bound or approved numerical contract establishes when that conversion is safe.

## Details

`fn compute_msg_to_child()` [`packages/treetime/src/partition/marginal_passes.rs#L185-L246`](../../packages/treetime/src/partition/marginal_passes.rs#L185-L246) divides and normalizes cavity profiles, then inserts them into `msg_to_child.variable` and decrements `fixed_counts`.

Let $\varepsilon$ denote a proposed probability threshold. A near-deterministic profile can retain nonzero probability outside its MAP state, so it is not identical to a fixed row. Demoting it under an $\varepsilon$ threshold would discard uncertainty and can change likelihood evaluation. The current explicit representation preserves the complete normalized profile.

## Investigation required

- Define a lossless equivalence criterion under the complete fixed-row representation, or derive an accumulated likelihood-error bound for approximate compression.
- Compare every affected sparse message and branch-likelihood contribution with the dense representation.
- Determine whether any approximate threshold can satisfy the project's numerical contract across repeated inference passes.
- Obtain explicit approval before introducing a lossy threshold.

No implementation ticket is ready until the representation criterion and numerical contract are decided.
