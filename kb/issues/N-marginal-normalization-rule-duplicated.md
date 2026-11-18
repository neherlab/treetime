# Marginal normalization rule is implemented twice

The rule that normalizes a probability vector and returns its log normalization exists in two functions, one for each sequence representation:

- `fn normalize_sites_inplace()` [packages/treetime/src/partition/marginal/shared/normalize.rs#L11-L25](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L11-L25): dense mode, every row of an `(L, K)` profile, returns the per-row log norms; `fn normalize_inplace()` [#L7-L9](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L7-L9) sums them. Callers: [packages/treetime/src/partition/marginal/shared/pass.rs#L243](../../packages/treetime/src/partition/marginal/shared/pass.rs#L243) (`normalize_sites_inplace`, root), `normalize_inplace` at [#L403](../../packages/treetime/src/partition/marginal/shared/pass.rs#L403), [#L410](../../packages/treetime/src/partition/marginal/shared/pass.rs#L410)
- `fn normalize_1d_inplace()` [packages/treetime/src/partition/marginal/sparse/message.rs#L255-L265](../../packages/treetime/src/partition/marginal/sparse/message.rs#L255-L265): sparse mode, one site or fixed-state vector, returns `weight` times its log norm. Callers: `fn compute_msg_to_child()` [packages/treetime/src/partition/marginal/sparse/forward.rs#L264](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L264) and [#L312](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L312)

Both apply the same rule to one vector: when its sum is positive and finite, divide by the sum and take the log of the sum; otherwise fill the vector uniformly with `1 / K` and report negative infinity. A change to the rule in one copy does not reach the other mode, and dense mode is the independent check of sparse mode, so the two copies must agree.

## The callers already diverge on degenerate vectors

`fn forward_log_lh_add_normalization()` [packages/treetime/src/partition/marginal/shared/normalize.rs#L41-L47](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L41-L47) skips a normalization of negative infinity. The two modes reach it with different granularity:

- Sparse `compute_msg_to_child` passes each site and each fixed state separately, so a degenerate vector drops only its own term
- Dense `msg_to_child` (`pass.rs#L402-L405`) passes the sum over all rows, so one degenerate row turns the sum into negative infinity and drops the log norms of every row of the message
- Dense node profiles (`pass.rs#L243`, `#L410`) add the sum directly, so one degenerate row makes the node log-likelihood negative infinity

> [!WARNING]
> **Needs review.** `pass.rs#L243` calls `normalize_sites_inplace` and keeps per-site root log norms (`root_log_norms`), so the claim above holds for `#L410` only unless the root path sums them later.


On input with a degenerate position, `msg_to_child.log_lh` therefore differs between dense and sparse mode. Which result is correct is the open question of [M-marginal-normalize-neg-infinity-masks-total.md](M-marginal-normalize-neg-infinity-masks-total.md).

## Fix direction

One function that normalizes one vector (an `ArrayViewMut1<f64>` with a weight), used by both modes: dense calls it per row and keeps its row-order sum, sparse calls it per site. The degenerate-vector semantics follow the decision of [M-marginal-normalize-neg-infinity-masks-total.md](M-marginal-normalize-neg-infinity-masks-total.md); until that decision, the merge keeps each mode's current result.

The merge keeps outputs byte-identical when the summation order of each vector stays the same. In ndarray 0.16.1, `sum_axis` along the axis of smallest stride sums each lane with `sum()` (`src/numeric/impl_numeric.rs`, `fn sum_axis`), which is the same eightfold unrolled sum (`src/numeric_util.rs`, `fn unrolled_fold`) that the sparse `sum()` of a contiguous vector uses. For a row-major `(L, K)` profile the per-row sums of the two functions are therefore identical.

> [!IMPORTANT]
> **Investigation required.** Confirm that every dense profile reaching `normalize_inplace` is row-major (`dis.is_standard_layout()` at the three callers). For a column-major profile, `sum_axis(Axis(1))` adds whole columns in sequence instead, which sums each row in a different order and would change outputs in the last bits for alphabets above 7 states. Validate a merge with byte-identical smoke outputs against `rust` in dense and sparse mode.

## Related issues

- [M-marginal-normalize-neg-infinity-masks-total.md](M-marginal-normalize-neg-infinity-masks-total.md): the degenerate-vector semantics

## Proposed implementation plan

1. Move the one-vector rule into `packages/treetime/src/partition/marginal/shared/normalize.rs` as one function on `ArrayViewMut1<f64>` with a weight; call it per row from `fn normalize_sites_inplace()` and per site from the sparse callers in `packages/treetime/src/partition/marginal/sparse/forward.rs#L264` and `#L312`; delete `fn normalize_1d_inplace()` from `packages/treetime/src/partition/marginal/sparse/message.rs`
2. Keep the degenerate-vector result of each mode until [M-marginal-normalize-neg-infinity-masks-total.md](M-marginal-normalize-neg-infinity-masks-total.md) is decided
3. Tests in `packages/treetime/src/partition/marginal/shared/__tests__/test_marginal_shared.rs`: a positive vector, an all-zero vector, and a vector with an infinite entry, each with weight 1 and weight 3 (oracle: hand calculation of the sum and `ln` of the sum); a row-major `(L, K)` profile normalized per row equals the per-vector results bit for bit (oracle: the one-vector function)
4. Validate with `just lint-rs`, `just test-rs`, and `just smoke` byte-identical against `rust` in dense and sparse mode, as the Investigation block above requires
