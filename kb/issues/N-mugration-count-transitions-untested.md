# Mugration count_transitions lacks direct unit test

## Needs investigation

The `count_transitions` function drives GTR parameter estimation for mugration. Its transition count matrix, dwell times, and root state filtering are exercised only through end-to-end golden master tests. A direct unit test with hand-computed expected values would catch matrix accumulation bugs, swapped parent-child axes, or wrong root-state filtering that produce silently wrong GTR parameters while stable argmax assignments mask the error.

## Details

`fn count_transitions()` at [packages/treetime/src/gtr/refinement.rs#L78-L110](../../packages/treetime/src/gtr/refinement.rs#L78-L110) accumulates:

- `nij`: expected transition count matrix from `get_branch_mutation_matrix` over all edges
- `Ti`: dwell times per state
- `root_state`: one-hot from root posterior argmax (filtered when near-uniform)

The function is private to `gtr::refinement`. It calls `get_branch_mutation_matrix` and `accumulate_mutation_counts` from `ancestral/gtr_inference_dense.rs`, which have their own tests. But the composition (edge iteration, branch length clamping, root filtering, diagonal zeroing) is untested in isolation.

## Proposed test

Construct a 3-leaf tree with known traits and branch lengths, attach traits, run a backward pass, then call `count_transitions()` and verify `nij`, `Ti`, and `root_state` against hand-computed values from the 2-state GTR transition equations.

## Current coverage

- Indirect: `test_gm_mugration_outputs` golden master tests (2 passing, 5 ignored)
- Indirect: `test_execute_mugration_*` integration tests in `test_run.rs`
