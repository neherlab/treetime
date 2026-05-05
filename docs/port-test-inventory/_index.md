# Test Inventory - TreeTime v1 (Rust)

## Summary

| Domain                                     | Page                                                                                                                      |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------- |
| [Ancestral Reconstruction](ancestral.md)   | Fitch parsimony, marginal ML (dense/sparse), equivalence, analytical, Python parity                                       |
| [Timetree Inference](timetree.md)          | Coalescent, belief propagation, convergence, confidence, relaxed clock, polytomies                                        |
| [GTR Models](gtr.md)                       | Model construction, matrix exponentiation, eigendecomposition, inference, site-specific, site-rate variation, JC distance |
| [Branch Optimization](optimization.md)     | Coefficient extraction, Newton-Raphson, grid search, dense/sparse equivalence, damping, indels, topology cleanup          |
| [Clock Inference](clock.md)                | Regression, date constraints, rerooting, best root, filtering, pipeline                                                   |
| [Mugration](mugration.md)                  | Discrete marginal, golden master, iterative refinement, confidence                                                        |
| [Distribution Operations](distribution.md) | Core, multiply, divide, convolve, negation, scalar, bounds, scaled, Gaussian, quantile                                    |
| [Representation](representation.md)        | Substitution composition, normalization, discrete states, topology cleanup, payloads                                      |
| [Supporting Crates](supporting.md)         | BitSet128, Alphabet, FASTA, Newick, sequence ops, graph, prune, CLI convert                                               |
| [Validation & Smoke](validation.md)        | CLI smoke tests, timetree/coalescent validation examples, metrics                                                         |

---

## Quick Reference

| Test Category            | Location                                                                                                                                              | Key Algorithms                                          |
| ------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------- |
| Fitch parsimony          | [`ancestral/__tests__/test_fitch.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_fitch.rs)                                         | Backward/forward passes, gap handling                   |
| Marginal ML dense        | [`ancestral/__tests__/test_marginal_dense.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_dense.rs)                       | Felsenstein pruning, sum-product                        |
| Marginal ML sparse       | [`ancestral/__tests__/test_marginal_sparse.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_sparse.rs)                     | Fitch compression + marginal                            |
| Dense/sparse equivalence | [`ancestral/__tests__/test_marginal_dense_sparse_*.rs`](../../packages/treetime/src/commands/ancestral/__tests__/)                                    | Cross-validation                                        |
| Python parity            | [`ancestral/__tests__/test_python_parity.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_python_parity.rs)                         | v0/v1 comparison                                        |
| Coalescent               | [`timetree/coalescent/__tests__/`](../../packages/treetime/src/commands/timetree/coalescent/__tests__/)                                               | Kingman, skyline, Tc optimization, total log-likelihood |
| Belief propagation       | [`timetree/inference/__tests__/`](../../packages/treetime/src/commands/timetree/inference/__tests__/)                                                 | Backward/forward time distribution                      |
| Relaxed clock            | [`timetree/optimization/__tests__/test_relaxed_clock.rs`](../../packages/treetime/src/commands/timetree/optimization/__tests__/test_relaxed_clock.rs) | Per-branch rate variation                               |
| GTR models               | [`gtr/__tests__/test_gm_gtr.rs`](../../packages/treetime/src/gtr/__tests__/test_gm_gtr.rs)                                                            | JC69, K80, HKY85, TN93, GTR                             |
| Matrix exponentiation    | [`gtr/__tests__/test_prop_gtr_expqt.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_expqt.rs)                                            | exp(Qt) properties                                      |
| GTR inference            | [`gtr/infer_gtr/__tests__/`](../../packages/treetime/src/gtr/infer_gtr/__tests__/)                                                                    | Dense/sparse mutation counts                            |
| Site-rate variation      | [`gtr/__tests__/test_prop_gtr_site_rates.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_site_rates.rs)                                  | Per-site rate heterogeneity                             |
| Newton-Raphson           | [`optimize/__tests__/test_newton_convergence/`](../../packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/)                     | Branch length optimization                              |
| Grid search              | [`optimize/__tests__/test_grid_search/`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/)                                   | Parameter sweeps                                        |
| Clock regression         | [`clock/__tests__/test_clock_regression.rs`](../../packages/treetime/src/commands/clock/__tests__/test_clock_regression.rs)                           | Root-to-tip analysis                                    |
| Rerooting                | [`clock/__tests__/test_reroot.rs`](../../packages/treetime/src/commands/clock/__tests__/test_reroot.rs)                                               | Optimal root placement                                  |
| Distributions            | [`treetime-distribution/src/__tests__/`](../../packages/treetime-distribution/src/__tests__/)                                                         | Multiply, divide, convolve                              |
| Smoke tests              | [`dev/run-smoke-tests`](../../dev/run-smoke-tests)                                                                                                    | End-to-end CLI validation                               |

---

## Test Type Definitions

| Type               | Description                       | Characteristics                       |
| ------------------ | --------------------------------- | ------------------------------------- |
| **Unit**           | Single function/module tests      | Fast, isolated, specific assertions   |
| **Property-based** | Random input testing (proptest)   | Random inputs, invariant verification |
| **Parameterized**  | Input/output case tables (rstest) | Named cases, aligned table format     |
| **Golden-master**  | v0 Python reference comparison    | Fixture files, tolerance thresholds   |
| **Integration**    | Multi-module workflows            | Real datasets, end-to-end pipelines   |

---

## Ignored Tests

Tests marked `#[ignore]` with reasons:

| Test                                                   | File                                                                                                                                                    | Reason                                                                             |
| ------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| `test_gm_mugration_outputs_v1_divergence`              | [`test_gm_mugration.rs`](../../packages/treetime/src/commands/mugration/__tests__/test_gm_mugration.rs)                                                 | v0 parity: intentional v1 improvements (pseudo-count pi, root-state filtering)     |
| `test_gm_mugration_confidence_outputs`                 | [`test_gm_mugration.rs`](../../packages/treetime/src/commands/mugration/__tests__/test_gm_mugration.rs)                                                 | v0 parity: confidence profiles exceed 1e-2 at multiple nodes (D1/D2 divergence)    |
| `test_prop_marginal_dense_sparse_gap_free_consistency` | [`test_marginal_dense_sparse_prop.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_dense_sparse_prop.rs)                     | Investigate dense-sparse marginal log-likelihood divergence on certain GTR configs |
| `test_discrete_gamma_rates_mean_one`                   | [`site_rate_variation.rs`](../../packages/treetime/src/gtr/site_rate_variation.rs)                                                                      | Flaky: discrete_gamma_rates fails for some alpha/K combinations                    |
| `test_prop_discrete_gamma_rates_mean_one`              | [`site_rate_variation.rs`](../../packages/treetime/src/gtr/site_rate_variation.rs)                                                                      | Flaky: discrete_gamma_rates fails for some alpha/K combinations                    |
| `test_runner_coalescent_completes`                     | [`test_runner_coalescent.rs`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_runner_coalescent.rs)               | Golden master datasets not yet passing                                             |
| `test_gm_runner_marginal_dense`                        | [`test_gm_runner_marginal_dense.rs`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs) | Golden master datasets not yet passing                                             |
| `test_gm_optimize`                                     | [`test_gm_optimize.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_gm_optimize.rs)                                                    | Per-branch divergence with v0 fixture exceeds 10% relative tolerance               |

---

## Loose Tolerances

Assertions using 1e-4 or looser:

| Tolerance | Test                                               | File                                                                                                                                                               | Rationale                                             |
| --------- | -------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------- |
| 1e-1      | `test_prop_marginal_sparse_log_lh_root_invariance` | [`test_marginal_root_invariance_prop.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_root_invariance_prop.rs)                          | Fitch ambiguity: compression pattern varies with root |
| 3e-1      | `test_gm_runner_poisson`                           | [`test_gm_runner_poisson.rs`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs)                          | Root node dominates max diff                          |
| 9e-1      | `test_gm_runner_marginal_dense`                    | [`test_gm_runner_marginal_dense.rs`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs)            | Disabled datasets                                     |
| 9e-1      | `test_gm_runner_marginal_sparse`                   | [`test_gm_runner_marginal_sparse.rs`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs)          | Disabled datasets                                     |
| 2e-2      | `test_gm_mugration_confidence_zika`                | [`test_gm_mugration.rs`](../../packages/treetime/src/commands/mugration/__tests__/test_gm_mugration.rs)                                                            | Mugration confidence comparison                       |
| 1e-2      | `test_convergence_iterations`                      | [`test_convergence_iterations.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_iterations.rs)                        | Cross-method LH comparison                            |
| 1e-2      | `test_gtr_numerical_edge_extreme_parameters`       | [`test_gtr_numerical_edge_extreme_parameters.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_extreme_parameters.rs) | Extreme parameter edge case                           |
| 1e-2      | `test_gm_gtr_site_specific_infer`                  | [`test_gm_gtr_site_specific.rs`](../../packages/treetime/src/gtr/__tests__/test_gm_gtr_site_specific.rs)                                                           | Circular inference distortion                         |
| 1e-2      | `test_infer_gtr_site_specific_*`                   | [`test_site_specific.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_site_specific.rs)                                                               | Circular inference distortion                         |
| 1e-3      | `test_optimize_method_*`                           | [`test_optimize_method.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_optimize_method.rs)                                                       | Cross-method LH agreement                             |
| 1e-3      | `test_gm_gtr_site_specific_infer` W ratios         | [`test_gm_gtr_site_specific.rs`](../../packages/treetime/src/gtr/__tests__/test_gm_gtr_site_specific.rs)                                                           | W ratio tolerance for inference                       |

---

## Commented-Out Test Cases

Datasets disabled via commented-out `#[case]` annotations:

| Test                             | Disabled Cases                                                              | Reason                                              |
| -------------------------------- | --------------------------------------------------------------------------- | --------------------------------------------------- |
| `test_gm_runner_poisson`         | dengue_20, lassa_L_20, mpox_clade_ii_20, rsv_a_20, tb_20, zika_20           | Missing internal node times, leaf dates not refined |
| `test_gm_runner_marginal_dense`  | dengue_20, ebola_20, lassa_L_20, mpox_clade_ii_20, rsv_a_20, tb_20, zika_20 | Same + ebola golden master node key mismatch        |
| `test_gm_runner_marginal_sparse` | dengue_20, lassa_L_20, mpox_clade_ii_20, rsv_a_20, tb_20, zika_20           | Same as poisson                                     |

---

## File Index by Crate

| Crate                   | Test Directories                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `treetime`              | `commands/ancestral/__tests__/`, `commands/clock/__tests__/`, `commands/mugration/__tests__/`, `commands/optimize/__tests__/`, `commands/prune/__tests__/`, `commands/timetree/**/__tests__/`, `gtr/__tests__/`, `alphabet/__tests__/`, `io/__tests__/`, `seq/__tests__/`, `graph/__tests__/`, `representation/__tests__/`, `representation/algo/topology_cleanup/__tests__/` + inline tests in `representation/`, `commands/optimize/`, `commands/timetree/coalescent/`, [`commands/timetree/run.rs`](../../packages/treetime/src/commands/timetree/run.rs) |
| `treetime-distribution` | `__tests__/`, `distribution_core/__tests__/`, `distribution_ops/__tests__/`, `distribution_scaled/__tests__/`                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| `treetime-analytical`   | `__tests__/`, inline tests                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `treetime-ops`          | `__tests__/`, inline tests                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `treetime-grid`         | `__tests__/`, inline tests                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `treetime-primitives`   | `__tests__/`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `treetime-utils`        | `datetime/__tests__/`, `fmt/__tests__/`, `array/__tests__/`, inline tests in `interval/`, `iterator/`, `array/`                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| `treetime-io`           | `__tests__/`, inline tests in [`parse_delimited.rs`](../../packages/treetime-io/src/parse_delimited.rs)                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `treetime-cli`          | `convert/__tests__/`, inline tests                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `treetime-validation`   | `testing/metrics/**` inline tests                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |

---

## Test Coverage Focus Areas

### Algorithmic Correctness

- Marginal reconstruction normalization (posteriors sum to 1.0)
- Root invariance (Felsenstein's pulley principle)
- Idempotency (single-pass convergence)
- Dense/sparse equivalence

### Numerical Stability

- Extreme branch lengths (1e-10 to 1000.0)
- Extreme equilibrium frequencies (pi values 0.0001 to 0.9997)
- High mutation rates (mu up to 100.0)
- Combined extreme parameters

### v0 Python Parity

- Ancestral sequence reconstruction
- Node time inference
- GTR parameter inference
- Coalescent contributions

### Topology Coverage

- Binary trees
- Polytomies (multifurcations)
- Caterpillar trees (linear)
- Deep trees (10+ levels)
- Star trees

---

## Running Tests

```bash
# All tests
./dev/docker/run ./dev/dev t

# Specific test pattern
./dev/docker/run ./dev/dev t test_marginal

# Smoke tests (end-to-end CLI)
./dev/docker/run ./dev/dev br && ./dev/docker/run ./dev/run-smoke-tests .build/docker/release/treetime

# Validation examples
./dev/docker/run ./dev/dev E timetree_validation
./dev/docker/run ./dev/dev E coalescent_validation
```

---

## Known Test Limitations

1. **Golden-master tolerance variations**: Some datasets require looser tolerances due to:
   - BLAS drift between NumPy and ndarray (scales with sequence length)
   - Fitch ambiguity resolution differences between dense/sparse
   - Zero-length branch handling differences

2. **Disabled datasets in runner tests**: Several datasets disabled due to zero-length branches:
   - `dengue_20`, `lassa_L_20`, `rsv_a_20`, `tb_20`, `zika_20`

3. **Dense marginal backward pass**: Crashes on some datasets (documented in validation example)

---

## Known Test Deficiencies

Tracked in `docs/port-known-issues/`:

- [Test coverage gaps](../port-known-issues/N-test-coverage-gaps.md): zero-test production functions across all major subsystems, 4 ignored golden-master tests, missing property tests for ClockSet and Fitch invariants
- [Test quality deficiencies](../port-known-issues/N-test-quality-deficiencies.md): missing pretty_assertions (46+ instances), loop-over-cases anti-pattern (21 files), circular tests, weak assertions, grossly loose tolerances (epsilon = 10.0, 1e0), helper placement violations

---

## References

- Test framework: [rstest](https://docs.rs/rstest) for parameterized tests
- Property testing: [proptest](https://docs.rs/proptest) for random input generation
- Floating-point comparison: [approx](https://docs.rs/approx) with ULP-based tolerances
- Assertions: [pretty_assertions](https://docs.rs/pretty_assertions) for clear diffs
