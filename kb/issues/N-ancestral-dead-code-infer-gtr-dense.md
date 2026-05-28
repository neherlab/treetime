# Dead production code in gtr_inference_dense.rs

`infer_gtr_dense` and its helper `get_mutation_counts_dense` in `gtr_inference_dense.rs` carry `#[allow(dead_code)]` annotations. `infer_gtr_dense` has no production callers. `get_mutation_counts_dense` is used only in tests.

The dense GTR inference functionality overlaps with mugration's `count_transitions_discrete`. The active sparse GTR inference path is `gtr_inference.rs`.

## Locations

- `packages/treetime/src/ancestral/gtr_inference_dense.rs:28-43:` `fn infer_gtr_dense()` -- `#[allow(dead_code)]`, no production callers
- `packages/treetime/src/ancestral/gtr_inference_dense.rs:145:` `fn get_mutation_counts_dense()` -- `#[allow(dead_code)]`, used only in tests

## Impact

Dead code increases maintenance surface. The `#[allow(dead_code)]` annotations suppress compiler warnings that would otherwise flag the unused code.

## Action

Evaluate whether the dense GTR inference path is needed for future work. If not, remove the dead code. If the test usage of `get_mutation_counts_dense` is valuable, consider whether those tests should use the sparse path instead.
