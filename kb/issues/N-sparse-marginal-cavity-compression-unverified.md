# Sparse cavity-profile compression lacks an error contract

> [!IMPORTANT]
> **Investigation required.** No analytical error bound or lossless equivalence criterion establishes when cavity-profile compression is safe. Collect the evidence listed under "Investigation required"; a lossy threshold then needs explicit approval.

## Summary

`fn compute_msg_to_child()` stores normalized cavity profiles explicitly even when their MAP state matches the fixed-row fallback and their probability is near $1$. Compressing those profiles approximately could improve representation density, but no analytical error bound or approved numerical contract establishes when that conversion is safe.

## Details

`fn compute_msg_to_child()` [`packages/treetime/src/partition/marginal/sparse/forward.rs#L219-L321`](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L219-L321) divides and normalizes cavity profiles, then inserts them into `msg_to_child.variable` and subtracts them from the fixed counts of their states.

Let $\varepsilon$ denote a proposed probability threshold. A near-deterministic profile can retain nonzero probability outside its MAP state, so it is not identical to a fixed row. Demoting it under an $\varepsilon$ threshold would discard uncertainty and can change likelihood evaluation. The current explicit representation preserves the complete normalized profile.

## Investigation required

- Define a lossless equivalence criterion under the complete fixed-row representation, or derive an accumulated likelihood-error bound for approximate compression.
- Compare every affected sparse message and branch-likelihood contribution with the dense representation.
- Determine whether any approximate threshold can satisfy the project's numerical contract across repeated inference passes.
- Obtain explicit approval before introducing a lossy threshold.

Implementation waits until the representation criterion and numerical contract are decided.

## Proposed implementation plan

The cited location is current (`fn compute_msg_to_child()` at `forward.rs#L219`). The plan stops at the approval of a threshold.

1. Measure first: count, on the example datasets, the cavity profiles whose MAP state equals the fixed-row state and whose MAP probability exceeds `1 - 1e-6`, `1 - 1e-9`, `1 - 1e-12`, and the share of `msg_to_child.variable` entries they make. Without a material share, close the issue
2. Lossless option: demote a profile only when it equals the fixed row in every component within the round-off of its normalization (a few ulps). This needs no approval and keeps dense and sparse equal within `1e-10`
3. Lossy option: derive the per-edge log-likelihood error of a demotion at threshold ε and its sum over sites and passes, and present it to the user for approval as a `kb/decisions/` entry before any code change

Tests: dense versus sparse log-likelihood and posteriors on small trees, equal within `1e-10` for the lossless option. Oracle: dense inference. Validation: `just test-rs`, `just lint-rs`.
