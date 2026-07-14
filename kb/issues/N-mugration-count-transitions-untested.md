# Mugration transition counting lacks a complete discrete hand oracle

## Needs investigation

The `count_transitions` function drives GTR parameter estimation for mugration. Shared contract tests already exercise parent/child orientation, edge accumulation, branch-length scaling, root composition, diagonal zeroing, and dense/sparse equality through `fn MarginalData::count_transitions()` [`packages/treetime/src/partition/marginal_core.rs#L483-L521`](../../packages/treetime/src/partition/marginal_core.rs#L483-L521). Missing coverage is narrower: direct discrete-partition delegation, a complete hand-computed two-state result, and explicit near-uniform root filtering.

## Details

`MarginalData::count_transitions()` at `packages/treetime/src/partition/marginal_core.rs` accumulates:

- `nij`: expected transition count matrix from `get_branch_mutation_matrix` over all edges
- `Ti`: dwell times per state
- `root_state`: argmax of root posterior per row (filtered when near-uniform)

The function is called via the `TransitionCounting` trait. Dense and discrete partitions delegate to the same implementation. Existing tests in [`packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs#L78-L119`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs#L78-L119) and [`packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs#L432-L476`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs#L432-L476) cover most shared behavior, but do not provide the missing discrete whole-result oracle.

## Proposed test

Construct a 3-leaf tree with known traits and branch lengths, attach traits, run a backward pass, then call `count_transitions()` and verify `nij`, `Ti`, and `root_state` against hand-computed values from the 2-state GTR transition equations.

## Current coverage

- Direct shared-contract tests for transition accumulation and dense/sparse consistency.
- Indirect `test_gm_mugration_outputs` golden masters and `test_execute_mugration_*` integration tests.

## Related tickets

- [kb/tickets/test-mugration-transition-counting.md](../tickets/test-mugration-transition-counting.md)
