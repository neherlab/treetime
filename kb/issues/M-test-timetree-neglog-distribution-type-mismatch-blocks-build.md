# Timetree lib tests do not compile after the negative-log distribution switch

Timetree time distributions now store negative-log ordinates, so the production payload holds `Arc<Distribution<NegLog>>` [`packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L62-L81`](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L62-L81). Several `treetime` test helpers still build a bare `Distribution`, which defaults to the `Plain` policy. The two types do not match, so the `treetime` lib test target fails to compile.

## Mechanism

`cargo clippy --all-targets -p treetime` fails with 10 `mismatched types` errors. Each error reports `Distribution<NegLog>` expected against `Distribution<Plain>` found. The test helpers construct a distribution without a policy parameter and pass it to a setter that now requires `NegLog`, for example [`packages/treetime/src/timetree/inference/__tests__/test_forward_pass.rs#L230`](../../packages/treetime/src/timetree/inference/__tests__/test_forward_pass.rs#L230) and [`packages/treetime/src/timetree/inference/__tests__/test_forward_pass.rs#L266`](../../packages/treetime/src/timetree/inference/__tests__/test_forward_pass.rs#L266). The affected files are `test_forward_pass.rs`, `test_backward_pass.rs`, `test_branch_length_likelihood.rs`, `test_date_constraints.rs`, and `test_confidence_extract.rs`.

## Impact

The build error blocks the whole `treetime` lib test target. The golden-master runners `test_gm_runner_marginal_dense`, `test_gm_runner_marginal_sparse`, and `test_gm_runner_poisson` cannot run. These runners are the acceptance criterion for the log-space distribution work, together with the large-dataset regression that must produce a finite positional log-likelihood and a date on every node. The workspace test command cannot pass until the helpers construct `Distribution<NegLog>`.

## Related

- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md): the validation plan names these golden masters as the acceptance criterion.
- [M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md): records the golden-master acceptance target that this build failure blocks.
