# Each edge contribution stores its own copy of the GTR model

`get_coefficients()` of both the dense and the sparse path stores `gtr.clone()` in every `PartitionContribution` it builds [packages/treetime/src/partition/optimize/dense.rs#L11](../../packages/treetime/src/partition/optimize/dense.rs#L11), [packages/treetime/src/partition/optimize/sparse.rs#L83](../../packages/treetime/src/partition/optimize/sparse.rs#L83). A contribution is built for every edge in every branch-length optimization round and every timetree iteration, so the model's matrices (`W`, `v`, `v_inv`) and vectors (`pi`, `eigvals`, optional `site_rates`) are copied once per edge each time.

## Impact

Negligible per edge compared with the per-site work, but it adds several heap allocations per edge per iteration, and every contribution of one partition holds an identical model.

## Proposed solution

Borrow the model instead of copying it: give `PartitionContribution` a lifetime and store `&GTR`, or share it through `Arc<GTR>` when a contribution must outlive the partition. The contribution reads the model only for `rate_eigvals` (computed at construction) and `certified_decay_rate()`.
