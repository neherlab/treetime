# Add pretty_assertions, trace attributes, and standardize approx wrappers

## Missing pretty_assertions (8 files, 46+ raw assert_eq!)

1. `commands/ancestral/__tests__/test_fitch.rs` (14)
2. `commands/clock/__tests__/test_date_constraints.rs` (14)
3. `commands/clock/__tests__/test_clock_filter.rs` (2)
4. `commands/clock/__tests__/test_clock_dengue100.rs` (1)
5. `commands/ancestral/__tests__/test_marginal_consistency.rs` (3)
6. `seq/__tests__/test_find_char_ranges.rs` (4)
7. `seq/__tests__/test_div.rs` (7)
8. `commands/mugration/__tests__/test_comment_output.rs` (1)

## Missing #[trace] on rstest (5 files, 17 blocks)

1. `gtr/__tests__/test_gm_gtr.rs` (7)
2. `gtr/__tests__/test_gm_gtr_site_specific.rs` (3)
3. `gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_parameterized.rs` (2)
4. `representation/__tests__/test_partition_marginal_sparse.rs` (1)
5. `seq/__tests__/test_find_char_ranges.rs` (4)

## Mixed approx wrappers (4 files)

Files using both `approx::assert_ulps_eq!` and `pretty_assert_ulps_eq!`:

1. `commands/ancestral/__tests__/test_marginal_dense.rs`
2. `commands/ancestral/__tests__/test_marginal_sparse.rs`
3. `commands/ancestral/__tests__/test_marginal_consistency.rs`
4. `commands/clock/__tests__/test_clock_regression.rs`

## Related issues

- Source: [N-test-quality-deficiencies.md](../issues/N-test-quality-deficiencies.md) -- delete after full resolution
