# Marginal normalization rule is implemented twice

The rule that normalizes a probability vector and returns its log normalization exists in two functions, one for each sequence representation:

- `fn normalize_inplace()` [packages/treetime/src/partition/marginal/shared/normalize.rs#L6-L21](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L6-L21): dense mode, every row of an `(L, K)` profile, returns the sum of the row log norms. Callers: [packages/treetime/src/partition/marginal/shared/pass.rs#L193](../../packages/treetime/src/partition/marginal/shared/pass.rs#L193), [#L302](../../packages/treetime/src/partition/marginal/shared/pass.rs#L302), [#L309](../../packages/treetime/src/partition/marginal/shared/pass.rs#L309)
- `fn normalize_1d_inplace()` [packages/treetime/src/partition/marginal/sparse/message.rs#L210-L219](../../packages/treetime/src/partition/marginal/sparse/message.rs#L210-L219): sparse mode, one site or fixed-state vector, returns `weight` times its log norm. Callers: `fn compute_msg_to_child()` [packages/treetime/src/partition/marginal/sparse/forward.rs#L207](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L207) and [#L220](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L220)

Both apply the same rule to one vector: when its sum is positive and finite, divide by the sum and take the log of the sum; otherwise fill the vector uniformly with `1 / K` and report negative infinity. A change to the rule in one copy does not reach the other mode, and dense mode is the independent check of sparse mode, so the two copies must agree.

## The callers already diverge on degenerate vectors

`fn forward_log_lh_add_normalization()` [packages/treetime/src/partition/marginal/shared/normalize.rs#L37-L43](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L37-L43) skips a normalization of negative infinity. The two modes reach it with different granularity:

- Sparse `compute_msg_to_child` passes each site and each fixed state separately, so a degenerate vector drops only its own term
- Dense `msg_to_child` (`pass.rs#L302-L305`) passes the sum over all rows, so one degenerate row turns the sum into negative infinity and drops the log norms of every row of the message
- Dense node profiles (`pass.rs#L193`, `#L309`) add the sum directly, so one degenerate row makes the node log-likelihood negative infinity

On input with a degenerate position, `msg_to_child.log_lh` therefore differs between dense and sparse mode. Which result is correct is the open question of [M-marginal-normalize-neg-infinity-masks-total.md](M-marginal-normalize-neg-infinity-masks-total.md).

## Fix direction

One function that normalizes one vector (an `ArrayViewMut1<f64>` with a weight), used by both modes: dense calls it per row and keeps its row-order sum, sparse calls it per site. The degenerate-vector semantics follow the decision of [M-marginal-normalize-neg-infinity-masks-total.md](M-marginal-normalize-neg-infinity-masks-total.md); until that decision, the merge keeps each mode's current result.

The merge keeps outputs byte-identical when the summation order of each vector stays the same. In ndarray 0.16.1, `sum_axis` along the axis of smallest stride sums each lane with `sum()` (`src/numeric/impl_numeric.rs`, `fn sum_axis`), which is the same eightfold unrolled sum (`src/numeric_util.rs`, `fn unrolled_fold`) that the sparse `sum()` of a contiguous vector uses. For a row-major `(L, K)` profile the per-row sums of the two functions are therefore identical.

> [!IMPORTANT]
> **Investigation required.** Confirm that every dense profile reaching `normalize_inplace` is row-major (`dis.is_standard_layout()` at the three callers). For a column-major profile, `sum_axis(Axis(1))` adds whole columns in sequence instead, which sums each row in a different order and would change outputs in the last bits for alphabets above 7 states. Validate a merge with byte-identical smoke outputs against `rust` in dense and sparse mode.

## Related issues

- [M-marginal-normalize-neg-infinity-masks-total.md](M-marginal-normalize-neg-infinity-masks-total.md): the degenerate-vector semantics; its source location `marginal_core.rs:175` is stale, the function is now in `partition/marginal/shared/normalize.rs`
- [N-ancestral-marginal-array-kernels-allocate.md](N-ancestral-marginal-array-kernels-allocate.md): allocations of the same marginal kernels
