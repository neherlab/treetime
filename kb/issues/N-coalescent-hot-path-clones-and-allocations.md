# Coalescent log weighting and model construction allocate avoidable buffers

The coalescent log-weighting operation allocates intermediate arrays once per internal node per backward pass, and the coalescent model clones invariant arrays once per timetree iteration.

## Log-weighting intermediate buffers

`fn distribution_multiply_by_fn()` [`packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L9-L44`](../../packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L9-L44) materializes the time grid with `distribution.t()` and collects a separate `weights` array before adding it to the ordinates. The weights could be added to the ordinates in place while walking the grid, without the time and weight arrays. The function runs once per internal node per backward pass and once per leaf with an active coalescent model.

## Model construction clones invariant arrays

`fn CoalescentModel::new()` [`packages/treetime/src/coalescent/coalescent.rs#L19-L26`](../../packages/treetime/src/coalescent/coalescent.rs#L19-L26) clones the lineage-count function and the Tc distribution. The model is constructed once per timetree iteration, so this is not per node, but the cloned values are invariant within the iteration and could be borrowed instead.

## Required behavior

Replace deep clones with borrows where ownership allows, and add the weights to the ordinates in place, keeping the arithmetic of every element unchanged so outputs stay byte-identical. Validate improvements with focused allocation benchmarks.
