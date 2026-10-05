# Commented-out code obscures test coverage and supported behavior

Maintained Rust sources keep disabled parameterized test cases and an unused `Seq::splice()` sketch as commented-out code. Each block starts with a `// TODO:` line that states why the code is disabled and, where one exists, the kb issue that blocks it:

- [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L35`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L35)
- [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L35`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L35)
- [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L26`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L26)
- [`packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L15`](../../packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L15) (two tests)
- [`packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs#L69`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs#L69)
- [`packages/treetime/src/gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs#L34`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs#L34)
- [`packages/treetime-primitives/src/seq.rs#L147`](../../packages/treetime-primitives/src/seq.rs#L147)

The `no_comments` lint keeps these blocks because their first line starts with the `TODO` marker (`dylint.toml`, section `[no_comments]`).

## Failure mode

- Commented `#[case]` rows are never collected, reported, or counted as ignored by the test runner, so test listings understate the datasets that lack coverage.
- Cases disabled as slow have no runner tier that executes them, so they rot silently when the fixtures or the code under test change.
- Commented production functions can silently become incompatible with the surrounding types while appearing to document intended capability.

## Open question

Whether a disabled case should move to an executable form (an `#[ignore]` case or a slow-test tier that CI runs), so that the runner reports it and compilation keeps it in sync with the code.
