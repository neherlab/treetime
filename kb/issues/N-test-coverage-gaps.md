# Test coverage gaps across production functions

## Summary

Systematic test coverage gaps across major subsystems: timetree inference, clock command, coalescent, ancestral, optimize, mugration, prune, representation, GTR, and foundation modules. Four golden-master tests are gated with `#[ignore]`.

## Ignored golden-master tests

- `timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs:40:` `#[ignore = "golden master datasets not yet passing"]`
- `timetree/inference/__tests__/test_gm_runner/test_runner_coalescent.rs:37:` `#[ignore = "golden master datasets not yet passing"]`
- `ancestral/__tests__/test_marginal_dense_sparse_prop.rs:75:` `#[ignore]` with `max_relative=1e-5`. Related: `M-ancestral-dense-sparse-divergence.md`
- `optimize/__tests__/test_gm_optimize.rs:70:` `#[ignore]`. Related: `M-optimize-gm-per-branch-divergence.md`

## Zero-test production functions

### Timetree inference

- `fn propagate_distributions_forward`: zero dedicated unit tests (forward pass tested only through GM)
- `fn create_poisson_branch_distributions`: no unit tests for edge cases
- `fn load_input_data` / `fn initialize_partitions`: no direct tests
- `fn run_refinement_iteration`: no tests
- `fn build_branch_distributions()`: untestable (`todo!()` body)
- `fn run_timetree_estimation()`: no branch coverage for input mode, confidence, skyline, rerooting, failure paths

### Clock command

- `fn run_clock()`: no end-to-end CLI test
- `fn clock_regression_forward`: no direct unit test
- `fn load_date_constraints()`: validation failure paths untested
- `fn write_clock_model()`, CSV writers, RTT chart writers: untested at serialized-output level

### Timetree optimization and output

- `fn apply_outlier_bad_branches` / `fn report_bad_branches`: zero tests
- `fn collect_outliers()` and `fn report_bad_branches()`: unverified
- `fn compute_rate_susceptibility()`: no integration coverage for upper/lower/restored rate passes
- `fn write_confidence_intervals()`: untested for TSV serialization
- `fn write_node_dates()`, `fn plot_root_to_tip()`, `fn plot_time_tree()`: unimplemented and untested

### Coalescent

- `fn collect_tree_events()` error paths untested (multiple roots, missing time distributions, non-finite present time)
- `fn collect_coalescent_edges()` and `fn sum_coalescent_cost()`: no edge-case coverage
- `fn compute_node_contributions()` and `fn compute_internal_contribution_single()`: no tests for `tc_dist.eval()` failure
- `fn compute_total_neg_log_lh()` and `fn optimize_skyline()`: no analytical or golden-master tests

### Ancestral command

- CLI entrypoint and `MethodAncestral` parsing: no end-to-end coverage
- Stdin FASTA path: no multi-record test
- `fn get_common_length()` error branches: untested
- `fn write_graph()` output files: untested for existence and parse-back

### Optimize command

- `fn apply_initial_guess_mode()`: no test for finite negative branch lengths in `Never` mode
- `fn run_optimize()`: no integration test proving negative branch lengths rejected
- `fn OptimizationContribution::from_sparse()` and `fn get_coefficients()`: not exercised through real sparse fixtures

### Mugration and prune

- `fn run_homoplasy()`: completely unverified (body is `unimplemented!()`)
- Mugration file-I/O wrappers: no integration coverage
- `fn optimize_gtr_rate()` and `fn refine_gtr_iterative()`: no tests for no-bracket path, backward-pass failure, rollback
- `fn run_prune`: no end-to-end test

### Representation module

- `fn combine_messages()`: zero direct unit tests (coverage indirect through integration only)
- `fn reconcile_topology()`: no tests (dense and sparse implementations)
- `fn fix_branch_length()`: no direct tests for clamp threshold, very short branches, `seq_length == 0`
- Sparse classification and reconstruction paths: no direct gap, unknown, ambiguity, parity, reroot-invariance tests
- `fn gather_points`: no test

### GTR

- `fn GTR::new()` invalid-input handling: no tests assert error returns vs panic
- Nucleotide constructors: no tests rejecting non-nucleotide alphabets
- `fn write_gtr_json()`: tests only check filename and existence, not JSON payload content
- `fn jtt92`: no direct regression coverage for 20-state empirical model

### Foundation

- `fn AlphabetConfig::validate()`: no direct test for `unknown` inside ambiguous value set
- `cli::rtt_chart` functions: no tests
- `timetree_validation.rs` functions: no tests for overlapping, disjoint, empty-overlap maps
- `seq::indel::InDel`: no dedicated constructor-boundary, inversion, formatting coverage

### Other

- No tests for `fn Sub::from_str`, `fn parse_pos`, validators at `seq/mutation.rs`
- No tests for `enum AlphabetName::AaNoStop` at `alphabet.rs`
- `fn count_sequence_changes`, `fn compute_rate_susceptibility`, `fn write_confidence_intervals`: zero direct tests
- `fn evaluate_site_contributions`: no direct unit test
- `trait BranchTopology` blanket impl: untested
- `fn propagate_raw`: not directly tested

## Missing property tests

### No property tests for ClockSet algebraic identities

`packages/treetime/src/payload/clock_set.rs:53-172:`

`+`, `-`, `+=`, `-=`, `fn propagate_averages` lack property test coverage for algebraic identities (associativity, commutativity, identity element).
