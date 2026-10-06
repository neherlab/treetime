# Test helpers defined outside the helpers module

Test files declare helper functions (fixtures, builders, assertion helpers) as top-level items of `mod tests`, mixed with or placed before the test functions. The project convention puts every helper after all tests, inside a nested `mod helpers { ... }` (property generators inside `mod generators { ... }`), with the tests importing them through `use self::helpers::{...}`.

## Impact

- Helpers and tests share one namespace, so a reader cannot tell a test from a helper without reading attributes
- New tests copy the pattern of the file they extend, so the deviation spreads to new code
- The `use self::helpers::{...}` line lists the shared setup that the tests depend on; without it, this dependency is implicit

## Example

[`packages/treetime-graph/src/__tests__/test_tree_view.rs`](../../packages/treetime-graph/src/__tests__/test_tree_view.rs) shows the required structure: `use self::helpers::{graph_from_parents, names};` at the top of `mod tests`, the tests, then `mod helpers` with its own explicit imports and `pub(super)` helpers.

## Locations

Files with at least one non-test function that is a direct item of `mod tests`, before any `mod helpers` or `mod generators`:

- `packages/app-commands/src/__tests__/test_examples_download.rs`
- `packages/app-commands/src/__tests__/test_support.rs`
- `packages/treetime-distribution/src/distribution_ops/__tests__/test_multiply.rs`
- `packages/treetime-distribution/src/distribution_ops/__tests__/test_prop_multiply_child_order.rs`
- `packages/treetime-distribution/src/distribution_ops/__tests__/test_prop_overlap.rs`
- `packages/treetime-grid/src/__tests__/hard_approach_law.rs`
- `packages/treetime-grid/src/__tests__/soft_tail_law.rs`
- `packages/treetime-utils/src/array/__tests__/test_batched.rs`
- `packages/treetime-utils/src/datetime/__tests__/parse_uncertain_date.rs`
- `packages/treetime/src/__tests__/test_error.rs`
- `packages/treetime/src/alphabet/__tests__/test_alphabet_config.rs`
- `packages/treetime/src/ancestral/__tests__/prop_marginal_support.rs`
- `packages/treetime/src/ancestral/__tests__/test_dense_completeness.rs`
- `packages/treetime/src/ancestral/__tests__/test_fitch_gap_sub_conflict.rs`
- `packages/treetime/src/ancestral/__tests__/test_fitch_indel.rs`
- `packages/treetime/src/ancestral/__tests__/test_fitch_sub.rs`
- `packages/treetime/src/ancestral/__tests__/test_marginal_analytical/test_marginal_analytical_support.rs`
- `packages/treetime/src/ancestral/__tests__/test_marginal_consistency.rs`
- `packages/treetime/src/ancestral/__tests__/test_marginal_dense.rs`
- `packages/treetime/src/ancestral/__tests__/test_marginal_dense_sparse_example.rs`
- `packages/treetime/src/ancestral/__tests__/test_marginal_idempotency_example.rs`
- `packages/treetime/src/ancestral/__tests__/test_marginal_map_deviation.rs`
- `packages/treetime/src/ancestral/__tests__/test_marginal_normalization_example.rs`
- `packages/treetime/src/ancestral/__tests__/test_marginal_normalization_prop.rs`
- `packages/treetime/src/ancestral/__tests__/test_marginal_sparse.rs`
- `packages/treetime/src/ancestral/__tests__/test_marginal_stability/test_marginal_stability_support.rs`
- `packages/treetime/src/ancestral/__tests__/test_marginal_tip_reconstruction.rs`
- `packages/treetime/src/ancestral/__tests__/test_mask.rs`
- `packages/treetime/src/ancestral/__tests__/test_python_parity.rs`
- `packages/treetime/src/clock/__tests__/test_clock_model.rs`
- `packages/treetime/src/clock/__tests__/test_clock_regression.rs`
- `packages/treetime/src/coalescent/__tests__/helpers.rs`
- `packages/treetime/src/coalescent/__tests__/test_coalescent_model.rs`
- `packages/treetime/src/coalescent/__tests__/test_events.rs`
- `packages/treetime/src/coalescent/__tests__/test_gm_coalescent.rs`
- `packages/treetime/src/coalescent/__tests__/test_gm_total_lh.rs`
- `packages/treetime/src/coalescent/__tests__/test_lineage_dynamics.rs`
- `packages/treetime/src/coalescent/__tests__/test_skyline.rs`
- `packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs`
- `packages/treetime/src/gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs`
- `packages/treetime/src/gtr/infer_gtr/__tests__/test_dense.rs`
- `packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs`
- `packages/treetime/src/mugration/__tests__/test_run.rs`
- `packages/treetime/src/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_invariants.rs`
- `packages/treetime/src/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_prop_invariants.rs`
- `packages/treetime/src/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_support.rs`
- `packages/treetime/src/optimize/__tests__/test_convergence/test_convergence_support.rs`
- `packages/treetime/src/optimize/__tests__/test_dense_edge_subs.rs`
- `packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_support.rs`
- `packages/treetime/src/optimize/__tests__/test_dispatch_zero_boundary.rs`
- `packages/treetime/src/optimize/__tests__/test_gm_optimize.rs`
- `packages/treetime/src/optimize/__tests__/test_grid_search/test_grid_search_support.rs`
- `packages/treetime/src/optimize/__tests__/test_initial_guess_formula.rs`
- `packages/treetime/src/optimize/__tests__/test_initial_guess_gaps.rs`
- `packages/treetime/src/optimize/__tests__/test_initial_guess_gtr_messages.rs`
- `packages/treetime/src/optimize/__tests__/test_is_zero_branch_optimal.rs`
- `packages/treetime/src/optimize/__tests__/test_newton_convergence/test_newton_convergence_support.rs`
- `packages/treetime/src/optimize/__tests__/test_optimize_indel.rs`
- `packages/treetime/src/optimize/__tests__/test_root_preservation.rs`
- `packages/treetime/src/optimize/topology/__tests__/test_collapse_edge.rs`
- `packages/treetime/src/optimize/topology/__tests__/test_merge_shared_mutations.rs`
- `packages/treetime/src/optimize/topology/__tests__/test_prop_merge_shared_mutations.rs`
- `packages/treetime/src/partition/marginal/dense/__tests__/test_marginal_dense.rs`
- `packages/treetime/src/partition/marginal/sparse/__tests__/test_partition_marginal_sparse.rs`
- `packages/treetime/src/partition/marginal/sparse/__tests__/test_sparse_transition_counting.rs`
- `packages/treetime/src/reroot/__tests__/test_cost_function.rs`
- `packages/treetime/src/reroot/__tests__/test_orchestrate.rs`
- `packages/treetime/src/seq/__tests__/test_composition.rs`
- `packages/treetime/src/seq/__tests__/test_gap_fill.rs`
- `packages/treetime/src/timetree/__tests__/test_coalescent_timescale.rs`
- `packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs`
- `packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs`
- `packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_support.rs`
- `packages/treetime/src/timetree/optimization/__tests__/test_reroot.rs`
- `packages/treetime/src/timetree/optimization/polytomy/__tests__/test_apply.rs`
- `packages/treetime/src/timetree/optimization/polytomy/__tests__/test_sweep.rs`
- `packages/util-newick/src/__tests__/test_nexus.rs`
- `packages/util-newick/src/__tests__/test_prop_roundtrip.rs`
- `packages/util-newick/src/__tests__/test_write.rs`
- `packages/util-phyloxml/src/__tests__/test_phyloxml_write.rs`

## Fix

Per file: move each helper into `mod helpers` after the last test, mark the helpers that tests call `pub(super)`, give `mod helpers` explicit imports, and import the helpers at the top of `mod tests` with `use self::helpers::{...}`. Replace glob imports of the module under test with named imports in the same change.

## Validation

- `just test-rs` runs the same number of tests per crate before and after the move
- `just lint-rs` reports no unused imports or unused helpers
