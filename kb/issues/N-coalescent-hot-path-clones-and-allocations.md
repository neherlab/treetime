# Coalescent log weighting allocates avoidable buffers

The coalescent log-weighting operation allocates intermediate arrays once per internal node per backward pass.

## Log-weighting intermediate buffers

`fn distribution_multiply_by_fn()` [`packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L9-L44`](../../packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L9-L44) materializes the time grid with `distribution.t()` and collects a separate `weights` array before adding it to the ordinates. The weights could be added to the ordinates in place while walking the grid, without the time and weight arrays. The function runs once per internal node per backward pass and once per leaf with an active coalescent model.

## Required behavior

Add the weights to the ordinates in place, keeping the arithmetic of every element unchanged so outputs stay byte-identical. Validate improvements with focused allocation benchmarks.
