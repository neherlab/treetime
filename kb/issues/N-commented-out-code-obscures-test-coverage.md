# Commented-out code obscures test coverage and supported behavior

Maintained Rust sources keep disabled parameterized test cases and an unused `Seq::splice()` sketch as commented-out code. Each block starts with a `// TODO:` line that states why the code is disabled and, where one exists, the kb issue that blocks it:

- [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L12`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L12)
- [`packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L15`](../../packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L15) and [`packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L95`](../../packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L95)
- [`packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs#L76`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs#L76)
- [`packages/treetime/src/gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs#L33`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs#L33)
- [`packages/treetime-primitives/src/seq.rs#L142`](../../packages/treetime-primitives/src/seq.rs#L142)

The `no_comments` lint keeps these blocks because their first line starts with the `TODO` marker (`dylint.toml`, section `[no_comments]`). The project coding conventions allow this form explicitly: disabled code is kept as a marker note under a `// TODO:` line that names the reason and the blocking kb issue. The two slow-dataset blocks in `infer_gtr` name no kb issue.

## Failure mode

- Commented `#[case]` rows are never collected, reported, or counted as ignored by the test runner, so test listings understate the datasets that lack coverage.
- Cases disabled as slow have no runner tier that executes them, so they rot silently when the fixtures or the code under test change.
- Commented production functions can silently become incompatible with the surrounding types while appearing to document intended capability.

## Open question

> [!IMPORTANT]
> **Decision required.** The project coding conventions allow disabled code under a `// TODO:` marker line, and this form hides the disabled cases from the test runner. Options:
>
> - **Keep the current form**: commented `#[case]` rows under `// TODO:` lines, as the conventions allow
> - **Move disabled cases to `#[ignore]`**: the runner lists and counts them, and compilation keeps them in sync with the code; the conventions change to require this form for disabled tests
> - **Move slow cases to a slow-test tier that CI runs**: slow datasets execute regularly; cases blocked by a kb issue stay `#[ignore]`
>
> The `Seq::splice()` sketch is production code with no caller; it is outside these options and can stay as a note or be deleted.
