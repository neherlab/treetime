# Golden-master runner cases are disabled for five datasets

The timetree golden-master runner tests compare v1 node times with v0 fixtures. Their cases for `dengue_20`, `lassa_L_20`, `mpox_clade_ii_20`, `rsv_a_20`, and `tb_20` are commented out under `// TODO:` lines in all three test files, so these datasets have no runner coverage:

- [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L26-L33](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L26-L33)
- [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L35-L43](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L35-L43)
- [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L35-L43](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L35-L43)

The fixtures still hold the expected values of these datasets.

## Current state

The original cause is gone. The cases were disabled because the inference produced times only for the leaves: the backward pass multiplied child messages whose finite grids were disjoint, and the product collapsed to `Empty`. The current code dates every node: a `timetree` run on `data/dengue/20` writes `numdate` for all 39 nodes.

The cases stay disabled because every remaining runner test is `#[ignore]`d for a gap to v0 that the grid resolution may cause ([M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md)): on `flu_h3n2_20`, the dense runner differs from v0 by up to 0.070 years and the Poisson runner by up to 0.35 years. The disabled cases may fail for further reasons once they run.

The `ebola_20` cases of the dense and sparse runners are disabled for a separate fixture defect ([N-timetree-dense-golden-master-node-mismatch.md](N-timetree-dense-golden-master-node-mismatch.md)), and the `zika_20` cases for the header defect of [M-timetree-date-header-hash.md](M-timetree-date-header-hash.md).

## Fix

- Re-enable the five dataset cases in all three runner tests after the runner tests pass at their target tolerance (`epsilon = 1e-6`)
- Investigate each case that then fails on its own, and file a separate issue for each new cause

## Validation

The 15 runner cases (5 datasets in 3 modes) run without `#[ignore]` and pass at the target tolerance.
