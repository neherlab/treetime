# Mugration count_transitions lacks direct unit test

## Needs investigation

The `count_transitions` function drives GTR parameter estimation for mugration. Its transition count matrix, dwell times, and root state filtering are exercised only through end-to-end golden master tests. A direct unit test with hand-computed expected values would catch matrix accumulation bugs, swapped parent-child axes, or wrong root-state filtering that produce silently wrong GTR parameters while stable argmax assignments mask the error.

## Details

`MarginalData::count_transitions()` at `packages/treetime/src/partition/marginal_core.rs` accumulates:

- `nij`: expected transition count matrix from `get_branch_mutation_matrix` over all edges
- `Ti`: dwell times per state
- `root_state`: argmax of root posterior per row (filtered when near-uniform)

The function is called via the `TransitionCounting` trait. Dense and discrete partitions delegate to `MarginalData::count_transitions`. The helper functions `get_branch_mutation_matrix` and `accumulate_mutation_counts` in `gtr/infer_gtr/common.rs` have their own tests. But the composition (edge iteration, branch length clamping, root filtering, diagonal zeroing) for discrete/mugration partitions is untested in isolation.

## Proposed test

Construct a 3-leaf tree with known traits and branch lengths, attach traits, run a backward pass, then call `count_transitions()` and verify `nij`, `Ti`, and `root_state` against hand-computed values from the 2-state GTR transition equations.

## Current coverage

- Indirect: `test_gm_mugration_outputs` golden master tests (2 passing, 5 ignored)
- Indirect: `test_execute_mugration_*` integration tests in `test_run.rs`
