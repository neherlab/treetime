# Enable 8 ignored golden-master and tolerance tests

## Ignored tests

### Blocked by grid resolution (`M-timetree-branch-grid-uniform-resolution.md`)

- [packages/treetime/src/commands/timetree/inference/\_\_tests\_\_/test_gm_runner/test_gm_runner_marginal_dense.rs](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs) `#[ignore]` at 1e-6. Max error 0.92 years.
- [packages/treetime/src/commands/timetree/inference/\_\_tests\_\_/test_gm_runner/test_gm_runner_poisson.rs](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs) `#[ignore]` at 1e-6. Max error 0.27 years.

### Blocked by algorithmic issues

- [packages/treetime/src/commands/ancestral/\_\_tests\_\_/test_marginal_root_invariance_prop.rs](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_root_invariance_prop.rs) sparse variant `#[ignore]` at 1e-6. Max diff ~1e-2. Related: `M-ancestral-sparse-root-invariance.md`.
- [packages/treetime/src/commands/ancestral/\_\_tests\_\_/test_marginal_dense_sparse_prop.rs](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_dense_sparse_prop.rs) `#[ignore]` with `max_relative=1e-5`. Related: `M-ancestral-dense-sparse-divergence.md`.

### Blocked by intentional v1 divergence

- [packages/treetime/src/commands/mugration/\_\_tests\_\_/test_gm_mugration.rs](../../packages/treetime/src/commands/mugration/__tests__/test_gm_mugration.rs) confidence zika `#[ignore]` at 1e-6. Max error 1.34e-2 from D1/D2 improvements. Related: `M-mugration-iterative-gtr.md`.

### Blocked by midpoint-rule discretization

- [packages/treetime/src/coalescent/\_\_tests\_\_/test_integration.rs](../../packages/treetime/src/coalescent/__tests__/test_integration.rs) `#[ignore]` at 1e-6. Midpoint-rule error 1.6e-4 on hyperbolic integrand. Related: `N-timetree-coalescent-integration-grossly-loose-tolerance.md`.

### Other

- [packages/treetime/src/commands/timetree/inference/\_\_tests\_\_/test_gm_runner/test_runner_coalescent.rs](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_runner_coalescent.rs) `#[ignore = "golden master datasets not yet passing"]`
- [packages/treetime/src/commands/optimize/\_\_tests\_\_/test_gm_optimize.rs](../../packages/treetime/src/commands/optimize/__tests__/test_gm_optimize.rs) `#[ignore]`. Related: `M-optimize-gm-per-branch-divergence.md`.

## Related issues

- Source: [N-test-coverage-gaps.md](../issues/N-test-coverage-gaps.md)
