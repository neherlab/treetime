# Evaluate and remove dead infer_gtr_dense code

`infer_gtr_dense` has no production callers and carries `#[allow(dead_code)]`.

## Current state

`ancestral/gtr_inference_dense.rs:28-43` defines `infer_gtr_dense()` with `#[allow(dead_code)]`. Its helper `get_mutation_counts_dense` (line 145) is also marked dead but used in tests. The active GTR inference path is `gtr_inference.rs` (sparse).

## Target state

If dense GTR inference is not needed: delete `infer_gtr_dense` and `get_mutation_counts_dense`, migrate any valuable test assertions to use the sparse inference path. If dense GTR inference is planned: remove `#[allow(dead_code)]` and wire it into production code.

## Implementation

1. Check whether any planned feature (e.g., dense-mode GTR inference, site-specific rates) requires `infer_gtr_dense`
2. If not needed: delete `gtr_inference_dense.rs` entirely, update `mod.rs`
3. If test coverage from `get_mutation_counts_dense` is valuable, port those test cases to use `gtr_inference.rs` equivalents
4. If needed: remove `#[allow(dead_code)]`, add production callers, add integration tests

## Related issues

Source: `kb/issues/N-ancestral-dead-code-infer-gtr-dense.md`
