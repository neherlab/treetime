# `range_contains()` scans every range for each position

`fn range_contains()` ([packages/treetime-utils/src/interval/range.rs#L9-L15](../../packages/treetime-utils/src/interval/range.rs#L9-L15)) tests a position against a list of ranges with a linear scan. The sparse marginal passes call it for every candidate site against the missing-data and non-character ranges of a node or an edge:

- [packages/treetime/src/partition/marginal/sparse/backward.rs#L192](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L192)
- [packages/treetime/src/partition/marginal/sparse/forward.rs#L151](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L151), [#L236](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L236), [#L351-L353](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L351-L353)
- [packages/treetime/src/partition/marginal/mixture/contribution.rs#L124](../../packages/treetime/src/partition/marginal/mixture/contribution.rs#L124)

The cost grows with the number of candidate sites times the number of ranges. On `data/mpox/clade-ii/100` (sparse, mostly missing data) with `--gamma-categories 4`, the scan takes 6.5% of the optimize run in a `perf` profile; runs without gamma call the same passes.

> [!IMPORTANT]
> **Investigation required.** A binary search or an `IntervalSet` lookup needs sorted, disjoint ranges. Check that every caller passes ranges in that form, then measure the run time of optimize on `mpox/clade-ii/100` and `/500` with and without gamma before and after the change.

## Validation

- The range lookup no longer shows in the profile of the sparse passes, and the outputs of optimize do not change
