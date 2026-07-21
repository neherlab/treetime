# Coalescent backward pass and log weighting allocate avoidable buffers on hot paths

The backward pass and coalescent log-weighting operation create intermediate allocations on per-node and per-optimizer-evaluation hot paths.

## Backward pass deep-clone

At [`packages/treetime/src/timetree/inference/backward_pass.rs#L98`](../../packages/treetime/src/timetree/inference/backward_pass.rs#L98), the completed internal-node distribution is cloned into a new `Arc`. The distribution was just computed and is not shared at this point; a move or in-place `Arc::new()` would avoid the deep copy.

## Log-weighting intermediate buffers

`fn distribution_apply_neg_log_weight()` [`packages/treetime-distribution/src/distribution_ops/log_cost.rs#L18-L58`](../../packages/treetime-distribution/src/distribution_ops/log_cost.rs#L18-L58) allocates a `weights` array via `collect()`, a `neg_log` array from the subtraction, and a `scaled` array from `mapv`. It also calls `distribution.t()` which may allocate a new `Array1` for Range variants. This function runs once per internal node per backward pass and once per leaf with an active coalescent model.

## Model construction clones invariant arrays

`fn CoalescentModel::new()` [`packages/treetime/src/coalescent/coalescent.rs#L47-L54`](../../packages/treetime/src/coalescent/coalescent.rs#L47-L54) clones the lineage-count array and the Tc distribution. The model is constructed once per timetree iteration, so this is not per-node, but the cloned arrays are invariant within the iteration and could be borrowed instead.

## Required behavior

Replace deep clones with moves or borrows where ownership allows. Reduce intermediate allocations in `distribution_apply_neg_log_weight` by combining operations or reusing buffers. Validate improvements with focused allocation benchmarks.
