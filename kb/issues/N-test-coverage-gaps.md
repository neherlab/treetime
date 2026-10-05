# Test coverage gaps across production functions

## Summary

Systematic test coverage gaps span timetree inference, clock, coalescent, ancestral reconstruction, optimize, mugration, prune, representation, GTR, and foundation modules. Relevant golden-master, analytical, and property tests remain gated with `#[ignore]`.

## Ignored golden-master tests

- Marginal dense golden master [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L40](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L40): `#[ignore = "golden master datasets not yet passing"]`
- Coalescent runner golden master [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_runner_coalescent.rs#L37](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_runner_coalescent.rs#L37): `#[ignore = "golden master datasets not yet passing"]`
- Dense/sparse property test [packages/treetime/src/ancestral/__tests__/test_marginal_dense_sparse_prop.rs#L12](../../packages/treetime/src/ancestral/__tests__/test_marginal_dense_sparse_prop.rs#L12): `#[ignore]` with `max_relative=1e-5`. Related: [M-ancestral-dense-sparse-divergence.md](M-ancestral-dense-sparse-divergence.md)
- Optimize golden master [packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L70](../../packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L70): `#[ignore]`. Related: [M-optimize-gm-per-branch-divergence.md](M-optimize-gm-per-branch-divergence.md)

## Zero-test production functions

### Timetree inference

- `fn propagate_distributions_forward`: zero dedicated unit tests (forward pass tested only through GM)
- `fn load_input_data` / `fn initialize_partitions`: no direct tests
- `fn run_timetree_estimation()`: no branch coverage for input mode, confidence, skyline, rerooting, failure paths
- `fn build_covariation_clock_params()` [packages/treetime/src/timetree/params.rs#L64-L94](../../packages/treetime/src/timetree/params.rs#L64-L94): no unit tests. It encodes the v0 covariation variance $\left(\ell + s^2/L\right)/L$ ([packages/legacy/treetime/treetime/clock_tree.py#L277-L285](../../packages/legacy/treetime/treetime/clock_tree.py#L277-L285)) as `variance_factor = 1/L`, `variance_offset = 0`, `variance_offset_leaf = s²/L²`, where $\ell$ is the branch length, $s$ the tip slack, and $L$ the sequence length. Tests in `packages/treetime/src/timetree/__tests__/test_params.rs`, with expected values from the v0 formula, not from running the function:
  - `covariation=false` returns `None`
  - `covariation=true, seq_len=Some(1000), tip_slack=None`: `variance_factor=1e-3`, `variance_offset=0.0`, `variance_offset_leaf = s²/1000²` for the default tip slack. The default value is undecided: [M-timetree-covariation-tip-slack-default-differs-from-v0.md](M-timetree-covariation-tip-slack-default-differs-from-v0.md)
  - `covariation=true, seq_len=Some(500), tip_slack=Some(5.0)`: exact values
  - `covariation=true, aln=Some(records)`: sequence length derived from the alignment
  - `covariation=true, seq_len=None, aln=None`: error

### Clock command

- `fn run_clock()`: no end-to-end CLI test
- `fn clock_regression_forward`: no direct unit test
- `fn load_date_constraints()`: validation failure paths untested
- Clock model JSON output, clock CSV output, RTT chart writers: untested at serialized-output level

### Timetree optimization and output

- `fn report_outliers()` and `fn collect_outlier_records()`: zero tests
- `fn compute_rate_susceptibility()`: no integration coverage for the inferences at the upper, lower and central rates
- `fn write_confidence_intervals()`: untested for TSV serialization

### Coalescent

- `fn collect_tree_events()` error paths untested (multiple roots, missing time distributions, non-finite present time)
- `fn collect_tree_events()`: error paths for multiple roots and missing time distributions remain untested
- `fn optimize_skyline()`: no independent end-to-end numerical oracle

### Ancestral command

- CLI entrypoint and `MethodAncestral` parsing: no end-to-end coverage
- Stdin FASTA path: no multi-record test
- `fn get_common_length()` error branches: untested
- `fn write_graph()` output files: untested for existence and parse-back

### Optimize command

- `fn apply_initial_guess_mode()`: no test for finite negative branch lengths in `Never` mode
- `fn run_optimize()`: no integration test proving negative branch lengths rejected
- `fn OptimizationContribution::from_sparse()` and `fn get_coefficients()`: not exercised through real sparse fixtures
- Damping guard of `fn run()` [packages/treetime/src/optimize/pipeline.rs#L42-L47](../../packages/treetime/src/optimize/pipeline.rs#L42-L47): no test that `damping >= 1.0` and `damping < 0.0` return the invalid-parameter error

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
- GTR JSON output: command tests only check filename and existence, not JSON payload content
- `fn jtt92`: no direct regression coverage for 20-state empirical model

### Foundation

- `fn AlphabetConfig::validate()`: no direct test for `unknown` inside ambiguous value set
- `cli::rtt_chart` SVG and PNG chart writers: no tests
- `timetree_validation.rs` functions: no tests for overlapping, disjoint, empty-overlap maps
- `seq::indel::InDel`: no dedicated constructor-boundary, inversion, formatting coverage

### Other

- No tests for `fn Sub::from_str`, `fn parse_pos`, validators at `seq/mutation.rs`
- No tests for `enum AlphabetName::AaNoStop` at `alphabet.rs`
- `fn count_sequence_changes`, `fn compute_rate_susceptibility`, `fn write_confidence_intervals`: zero direct tests
- `fn evaluate_site_contributions`: no direct unit test
- `trait BranchTopology` blanket impl: untested
- `fn propagate_raw`: not directly tested

## Cross-cutting scientific coverage gaps

- Scientific oracle and property tests remain ignored across distribution, rerooting, coalescent, and inference paths.
- Graph-backed semantic output tests cover `timetree` more thoroughly than ancestral, clock, mugration, optimize, and prune.
- Parallel coverage concentrates on one sparse success case and does not establish error atomicity for marginal, optimize, or timetree passes.
- Skyline tests recompute the reported formula instead of invoking the optimizer objective.
- Coalescent initialization tests establish enum selection without checking topology-change sequencing or event completeness.
- Sentinel arithmetic is tested as isolated helpers without a multiple-impossible-factor cavity case.
- Fitch properties cover gap-free sequence reversal but not arbitrary column permutations, ambiguity, or exhaustive multifurcation scores.
- Benchmark/report tooling lacks an automated fixture for revision dimensions and requested worker counts.

Production defects remain in their domain issues; this issue owns the cross-cutting test matrix and re-enabling valid ignored tests.

## Missing property tests

### No property tests for ClockSet algebraic identities

`struct ClockSet` algebra and propagation [packages/treetime/src/clock/clock_set.rs#L55-L182](../../packages/treetime/src/clock/clock_set.rs#L55-L182)

Addition and `+=` lack coverage for associativity, commutativity, and the zero identity. Subtraction and `-=` require inverse-law properties such as `(a - b) + b = a`; subtraction is neither associative nor commutative. `fn propagate_averages` needs separately derived invariants for its valid input domain.

### No property tests for Fitch parsimony invariants

Score invariance under rerooting remains relevant. State-set containment must be conditioned on the Fitch recurrence: an intersection result is a subset of every child set, while a union result contains each disjoint child set. A single unconditional parent/child subset property is false.

## Potential solutions

- O1. Add tests in focused changes, one production ownership boundary at a time, after the corresponding behavior and oracle are defined.
- O2. Add coverage for every listed function and ignored suite in one change. This obscures distinct oracles and makes blocked production defects appear test-ready.

## Recommendation

Use O1. Keep this file as the coverage inventory, remove each entry when its tests land, and enable an ignored test only after its production or parity blocker is resolved.

## Readiness

The inventory as a whole is not ready for implementation. Focused property and domain tests require valid properties, explicit generators, and traceable oracles; ignored golden masters and unrelated zero-test functions require separate issues or resolved blockers.
