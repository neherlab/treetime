# Zero-length branches cause panic

Zero-length or near-zero branches in the input tree propagate through belief
propagation, producing distribution grids with non-uniform spacing. The
interpolation layer panics with `x array must be uniformly spaced`. Affected
datasets: dengue_20, lassa_L_20, mpox_clade_ii_20, rsv_a_20, tb_20. These
datasets are commented out in the timetree runner tests:

- [test_gm_runner_poisson.rs#L19](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L19)
- [test_gm_runner_marginal_dense.rs#L31](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L31)
- [test_gm_runner_marginal_sparse.rs#L32](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L32)

Related to [zero branch length clamping](N-core-branch-length-clamping.md) but
occurs at a different point in the pipeline - the clamping in
`fix_branch_length` applies during ancestral reconstruction, not during
timetree distribution construction.
