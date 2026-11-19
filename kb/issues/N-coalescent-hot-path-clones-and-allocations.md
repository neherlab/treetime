# Coalescent log weighting allocates avoidable buffers

The coalescent log-weighting operation allocates intermediate arrays on every call: a time array, a clone of the ordinates, and a separate weight array.

## Log-weighting intermediate buffers

`fn distribution_multiply_by_fn()` [`packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L8-L45`](../../packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L8-L45) materializes the time grid with `distribution.t()` (an owned `Array1`, L25), clones the ordinates of a `Distribution::Function` (L27), and collects a separate `weights` array (L33) before it adds the weights to the cloned ordinates (L34; the owned addition reuses the clone's buffer). The weights could be added to the ordinates in place while walking the grid, without the time and weight arrays. The callers are in `packages/treetime/src/timetree/inference/coalescent_factor.rs`; the function runs once per internal node per backward pass and once per leaf with an active coalescent model.

## Required behavior

Add the weights to the ordinates in place, keeping the arithmetic of every element unchanged so outputs stay byte-identical. Validate improvements with focused allocation benchmarks.

## Proposed implementation plan

1. In `fn distribution_multiply_by_fn()`, take the ordinates as before (one owned copy is needed, because the input is borrowed) and add `weight(t)` to each element in place with `ndarray::Zip` over the ordinates and the grid abscissas of the `Function` (read through the grid view of `f.grid_fn()`, without `distribution.t()`), propagating the first `Err` of `weight`
2. Keep `distribution.t()` only in the `Point` and `Range` branch, where the result is built from the times
3. Keep the element arithmetic `y + w` in the same order, so the outputs stay byte-identical
4. Tests: the existing cases of `packages/treetime-distribution/src/distribution_ops/__tests__/test_multiply_by_fn.rs` stay unchanged and must pass with exact equality (oracle: the outputs of the current implementation, which the change must not alter, plus the existing analytical expectations of that file)
5. Validation: `just check`, `just test-rs -E 'package(treetime-distribution)'`, then `just test-rs`; measure allocations and time of the backward pass with `dev/profile` on `data/sc2/4500` with `data/sc2/4500/timetree.yaml` (skyline coalescent) before and after the change, and report both
