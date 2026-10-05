# Mugration transition counting lacks a complete discrete hand oracle

## Needs investigation

Transition counting drives GTR parameter estimation for mugration. Shared contract tests already exercise parent/child orientation, edge accumulation, branch-length scaling, root composition, diagonal zeroing, and dense/sparse equality through `fn count_transitions_dense()` [`packages/treetime/src/partition/marginal/shared/data.rs#L16-L57`](../../packages/treetime/src/partition/marginal/shared/data.rs#L16-L57). Missing coverage is narrower: direct discrete-partition delegation, a complete hand-computed two-state result, and explicit near-uniform root filtering.

## Details

`fn count_transitions_dense()` at `packages/treetime/src/partition/marginal/shared/data.rs` accumulates:

- `nij`: expected transition count matrix from `get_branch_mutation_matrix` over all edges
- `Ti`: dwell times per state
- `root_state`: argmax of root posterior per row (skipped when near-uniform and `filter_uninformative_root` is set)

The function is called through `fn count_transitions()` of the `MarginalPasses` trait [`packages/treetime/src/partition/marginal/shared/update.rs#L39`](../../packages/treetime/src/partition/marginal/shared/update.rs#L39). Dense and discrete partitions delegate to the same implementation ([`packages/treetime/src/partition/marginal/dense/partition.rs#L280`](../../packages/treetime/src/partition/marginal/dense/partition.rs#L280), [`packages/treetime/src/partition/marginal/discrete/partition.rs#L136-L146`](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L136-L146)). Existing tests in [`packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs#L74-L119`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs#L74-L119) and [`packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs#L413-L461`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs#L413-L461) cover most shared behavior, but do not provide the missing discrete whole-result oracle.

## Proposed test

Construct a 3-leaf tree with known traits and branch lengths, attach traits, run a backward pass, then call `count_transitions()` through the discrete partition and verify `nij`, `Ti`, and `root_state` against hand-computed values from the 2-state GTR transition equations.

## Acceptance criteria

- A hand-computed two-state fixture verifies the full expected transition matrix, dwell-time vector, and root state
- The hand oracle distinguishes near-uniform root filtering and verifies discrete delegation; existing tests remain the oracle for parent/child orientation, diagonal zeroing, branch-length weighting, and dense/sparse equality
- The oracle does not call the production accumulation helpers under test

## Current coverage

- Direct shared-contract tests for transition accumulation and dense/sparse consistency.
- Indirect `test_gm_mugration_outputs` golden masters and `test_run_mugration_*` integration tests in `packages/treetime/src/mugration/__tests__/test_run.rs`.
