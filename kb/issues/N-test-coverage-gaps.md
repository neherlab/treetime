# Test coverage gaps across production functions

## Summary

Systematic test coverage gaps span timetree inference, clock, coalescent, ancestral reconstruction, optimize, mugration, prune, representation, GTR, and foundation modules. Golden-master, property, and dense/sparse consistency tests remain gated with `#[ignore]`.

## Ignored tests

Golden masters:

- Optimize per-branch golden master [packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L50](../../packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L50): `#[ignore]`. Related: [M-optimize-gm-per-branch-divergence.md](M-optimize-gm-per-branch-divergence.md)

Property, consistency, and command tests:

- Dense/sparse property test [packages/treetime/src/ancestral/__tests__/test_marginal_dense_sparse_prop.rs#L50](../../packages/treetime/src/ancestral/__tests__/test_marginal_dense_sparse_prop.rs#L50): `#[ignore]` with `max_relative = 1e-5`. Related: [M-ancestral-dense-sparse-divergence.md](M-ancestral-dense-sparse-divergence.md)
- Sparse root-invariance property test [packages/treetime/src/ancestral/__tests__/test_marginal_root_invariance_prop.rs#L27](../../packages/treetime/src/ancestral/__tests__/test_marginal_root_invariance_prop.rs#L27): `#[ignore]`. Related: [M-ancestral-sparse-root-invariance.md](M-ancestral-sparse-root-invariance.md)
- Sparse GTR refinement against dense [packages/treetime/src/gtr/__tests__/test_refinement.rs#L202](../../packages/treetime/src/gtr/__tests__/test_refinement.rs#L202): `#[ignore]`. Related: [M-gtr-sparse-transition-counts-diverge-from-dense.md](M-gtr-sparse-transition-counts-diverge-from-dense.md)
- Homoplasy dense/sparse comparison [packages/app-commands/src/commands/homoplasy/__tests__/test_dense_sparse.rs#L17](../../packages/app-commands/src/commands/homoplasy/__tests__/test_dense_sparse.rs#L17): `#[ignore]`. Related: [H-ancestral-sparse-dense-internal-residues-diverge.md](H-ancestral-sparse-dense-internal-residues-diverge.md)
- Ancestral ambiguity report [packages/app-commands/src/commands/ancestral/__tests__/test_report_ambiguous.rs#L86](../../packages/app-commands/src/commands/ancestral/__tests__/test_report_ambiguous.rs#L86): `#[ignore]`. Related: [M-ancestral-sparse-imputation-skips-tips-below-unknown-node.md](M-ancestral-sparse-imputation-skips-tips-below-unknown-node.md)

These ignored tests are the independent oracles of the project: the v0 golden masters and the comparisons of sparse with dense results. While they are ignored, a change that widens one of their gaps passes CI unnoticed.

> [!IMPORTANT]
> **Decision required.** How to guard the ignored reference comparisons until their gaps are fixed. The options are:
>
> - Keep them ignored, and rerun them by hand with `--run-ignored only` when a related change lands
> - Add a regression guard beside each one that asserts the currently measured difference as an upper bound, so a gap can shrink but not grow unnoticed. The guard is not a correctness test and goes away when the issue that owns the gap is fixed; the ignored test keeps the target tolerance

## Zero-test production functions

### Timetree inference

- `fn propagate_distributions_forward`: zero dedicated unit tests (forward pass tested only through GM)
- `fn load_input_data`: no direct tests
- `fn run_timetree_estimation()`: no branch coverage for input mode, confidence, skyline, rerooting, failure paths

### Clock command

- `fn run_clock()`: no end-to-end CLI test
- `fn clock_regression_forward`: no direct unit test
- `fn load_date_constraints()`: validation failure paths untested
- Clock model JSON output, clock CSV output: untested at serialized-output level
- SVG and PNG root-to-tip chart writers [packages/app-commands/src/rtt_chart.rs#L21-L52](../../packages/app-commands/src/rtt_chart.rs#L21-L52): no tests of the written files; only `fn gather_points()` has tests

### Timetree optimization and output

- `fn compute_rate_susceptibility()`: no integration coverage for the inferences at the upper, lower and central rates

### Coalescent

- `fn optimize_skyline()`: no independent end-to-end numerical oracle

### Ancestral command

- CLI entrypoint and `MethodAncestral` parsing: no end-to-end coverage
- Stdin FASTA path: no multi-record test
- `fn common_length()` error branches: untested

### Optimize command

- `fn apply_initial_guess_mode()`: no test for finite negative branch lengths in `Never` mode
- `fn run_optimize()`: no integration test proving negative branch lengths rejected
- `fn OptimizationContribution::from_sparse()` and `fn get_coefficients()`: not exercised through real sparse fixtures
- Damping guard of `fn run()` [packages/treetime/src/optimize/pipeline.rs#L43-L48](../../packages/treetime/src/optimize/pipeline.rs#L43-L48): no test that `damping >= 1.0` and `damping < 0.0` return the invalid-parameter error

### Mugration and prune

- Mugration file-I/O wrappers: no integration coverage
- `fn optimize_gtr_rate()`: no tests for no-bracket path, backward-pass failure, rollback
- `fn run_prune`: no end-to-end test

### Representation module

- `fn combine_messages()`: zero direct unit tests (coverage indirect through integration only)
- `fn reconcile_topology()`: no tests (dense and sparse implementations)
- `fn fix_branch_length()`: no direct tests for clamp threshold, very short branches, `seq_length == 0`
- Sparse classification and reconstruction paths: no direct gap, unknown, ambiguity, parity, reroot-invariance tests

### GTR

- `fn GTR::new()` invalid-input handling: no tests assert error returns vs panic
- Nucleotide constructors: no tests rejecting non-nucleotide alphabets
- GTR JSON output: command tests only check filename and existence, not JSON payload content
- `fn jtt92`: no direct regression coverage for 20-state empirical model

### Foundation

- `fn AlphabetConfig::validate()`: no direct test for `unknown` inside ambiguous value set
- `timetree_validation.rs` functions: no tests for overlapping, disjoint, empty-overlap maps
- `seq::indel::InDel`: no dedicated constructor-boundary, inversion, formatting coverage

### Other

- No tests for `fn parse_pos` and the validators in `seq/mutation.rs`
- `fn count_sequence_changes`, `fn compute_rate_susceptibility`: zero direct tests
- `fn evaluate_site_contributions`: no direct unit test
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

`struct ClockSet` algebra and propagation [packages/treetime/src/clock/clock_set.rs](../../packages/treetime/src/clock/clock_set.rs)

Addition and `+=` lack coverage for associativity, commutativity, and the zero identity. Subtraction and `-=` require inverse-law properties such as `(a - b) + b = a`; subtraction is neither associative nor commutative. `fn ClockSet::shift` should compose additively (`shift(a).shift(b) = shift(a + b)`) and commute with addition.

### No property tests for Fitch parsimony invariants

Score invariance under rerooting remains relevant. State-set containment must be conditioned on the Fitch recurrence: an intersection result is a subset of every child set, while a union result contains each disjoint child set. A single unconditional parent/child subset property is false.

## Potential solutions

- O1. Add tests in focused changes, one production ownership boundary at a time, after the corresponding behavior and oracle are defined.
- O2. Add coverage for every listed function and ignored suite in one change. This obscures distinct oracles and makes blocked production defects appear test-ready.

## Recommendation

Use O1. Keep this file as the coverage inventory, remove each entry when its tests land, and enable an ignored test only after its production or parity blocker is resolved.

## Readiness

The inventory as a whole is not ready for implementation. Focused property and domain tests require valid properties, explicit generators, and traceable oracles; ignored golden masters and unrelated zero-test functions require separate issues or resolved blockers.
