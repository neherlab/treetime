# `range_contains()` scans every range for each position

`fn range_contains()` ([packages/treetime-utils/src/interval/range.rs#L9-L15](../../packages/treetime-utils/src/interval/range.rs#L9-L15)) tests a position against a list of ranges with a linear scan. The sparse marginal passes call it for every candidate site against the missing-data and non-character ranges of a node or an edge:

- [packages/treetime/src/partition/marginal/sparse/backward.rs#L192](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L192)
- [packages/treetime/src/partition/marginal/sparse/forward.rs#L151](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L151), [#L236](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L236), [#L351-L353](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L351-L353)
- [packages/treetime/src/partition/marginal/mixture/contribution.rs#L124](../../packages/treetime/src/partition/marginal/mixture/contribution.rs#L124)
- Further callers in the sparse backward pass, the Fitch substitutions, the dense marginal and the topology collapse: [packages/treetime/src/partition/marginal/sparse/backward.rs#L193](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L193), [#L232](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L232), [packages/treetime/src/partition/fitch/sub.rs#L29](../../packages/treetime/src/partition/fitch/sub.rs#L29) (and L36, L138, L153, L160), [packages/treetime/src/partition/marginal/dense/partition.rs#L100](../../packages/treetime/src/partition/marginal/dense/partition.rs#L100) (and L218), [packages/treetime/src/optimize/topology/collapse.rs#L52](../../packages/treetime/src/optimize/topology/collapse.rs#L52)

The cost grows with the number of candidate sites times the number of ranges. On `data/mpox/clade-ii/100` (sparse, mostly missing data) with `--gamma-categories 4`, the scan takes 6.5% of the optimize run in a `perf` profile; runs without gamma call the same passes.

> [!IMPORTANT]
> **Investigation required.** A binary search or an `IntervalSet` lookup needs sorted, disjoint ranges. Check that every caller passes ranges in that form, then measure the run time of optimize on `mpox/clade-ii/100` and `/500` with and without gamma before and after the change.

## Validation

- The range lookup no longer shows in the profile of the sparse passes, and the outputs of optimize do not change

## Proposed implementation plan

1. Check the producers of `gaps`, `non_char` and `child_no_data` (the Fitch compression and the sparse backward pass) for sorted, disjoint, half-open ranges; if a producer does not guarantee it, normalize there (the owner of the data), not in the lookup
2. Replace the scan of `fn range_contains()` in `packages/treetime-utils/src/interval/range.rs` with a binary search (`partition_point` on the range starts, then a check of the end) on the sorted form, and keep `range_contains_iter()` only if a caller still needs an unsorted iterator
3. Tests in `packages/treetime-utils/src/interval/__tests__/`: an `rstest` table with positions before, at the start, inside, at the end (excluded) and after each range, and between ranges; a property test that the binary search matches the linear scan on random sorted disjoint ranges (oracle: the linear scan, an independent formula of membership)
4. Validate: `just test-rs`; optimize outputs on `mpox/clade-ii/100` and `/500` with and without `--gamma-categories 4` are byte-identical before and after (`just smoke` cases or `just r`), and the host profile (`dev/profile`) no longer shows the lookup
