# Enable 4 ignored golden-master tests

## Ignored tests

- `timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs:40:` `#[ignore = "golden master datasets not yet passing"]`
- `timetree/inference/__tests__/test_gm_runner/test_runner_coalescent.rs:37:` `#[ignore = "golden master datasets not yet passing"]`
- `ancestral/__tests__/test_marginal_dense_sparse_prop.rs:75:` `#[ignore]` with `max_relative=1e-5`. Related: `M-ancestral-dense-sparse-divergence.md`
- `optimize/__tests__/test_gm_optimize.rs:70:` `#[ignore]`. Related: `M-optimize-gm-per-branch-divergence.md`

## Related issues

- Source: [N-test-coverage-gaps.md](../issues/N-test-coverage-gaps.md) -- delete after full resolution
