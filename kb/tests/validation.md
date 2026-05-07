# Validation & Smoke Tests

[Back to index](README.md)

## Summary

| Category                                                | Type          |
| ------------------------------------------------------- | ------------- |
| [Smoke tests (CLI)](#smoke-tests)                       | Integration   |
| [Timetree validation](#timetree-validation)             | Golden-master |
| [Coalescent validation](#coalescent-validation)         | Golden-master |
| [treetime-validation crate](#treetime-validation-crate) | Unit          |

---

## Smoke Tests

**Script:** [`dev/run-smoke-tests`](../../dev/run-smoke-tests)

End-to-end CLI validation running the treetime binary across datasets and flag combinations. Executed in parallel via Python multiprocessing.

### Breakdown

| Command   | Modes                                      |
| --------- | ------------------------------------------ |
| ancestral | parsimony, marginal-sparse, marginal-dense |
| clock     | default                                    |
| optimize  | sparse (dense disabled)                    |
| prune     | default                                    |
| timetree  | varied flag combinations                   |
| help      | --help per subcommand                      |

### Datasets Used

`dengue/{20,100}`, `ebola/{20,100,362}`, `flu/h3n2/{20,200}`, `lassa/L/{20,50,200}`, `rsv/a/{20,100}`, `tb/{20,100,149}`, `zika/{20,86}`

Additional timetree-only: `mpox/clade-ii/20`

### Timetree Flag Coverage

| Category           | Flags tested                                                    |
| ------------------ | --------------------------------------------------------------- |
| Basic              | defaults                                                        |
| Clock rate         | `--clock-rate=0.003`                                            |
| Branch length mode | `--branch-length-mode=input`                                    |
| Time marginal      | `--time-marginal={always,only-final}`, `--confidence`           |
| Polytomies         | `--keep-polytomies`, `--resolve-polytomies`                     |
| Coalescent         | `--coalescent=10.0`, `--coalescent-opt`, `--coalescent-skyline` |
| Reroot             | `--reroot={min-dev,oldest,least-squares}`, `--keep-root`        |
| Clock filter       | `--clock-filter={0,5.0}`                                        |
| Covariation        | `--covariation`                                                 |
| Ancestral method   | `--method-anc=parsimony`                                        |
| GTR model          | `--gtr=jc69`                                                    |
| Max iterations     | `--max-iter={0,1,4}`                                            |
| Tip slack          | `--tip-slack=0.1`                                               |
| Combinations       | multi-flag combinations                                         |
| Known failures     | `--vary-rate`, `--relax`, `--dense=true`, crash-prone datasets  |

---

## Timetree Validation

**Test:** [`packages/treetime/examples/timetree_validation.rs`](../../packages/treetime/examples/timetree_validation.rs)

**Impl:**

- [`packages/treetime/src/commands/timetree/inference/runner.rs`](../../packages/treetime/src/commands/timetree/inference/runner.rs)
- [`packages/treetime/src/commands/timetree/inference/backward_pass.rs`](../../packages/treetime/src/commands/timetree/inference/backward_pass.rs)
- [`packages/treetime/src/commands/timetree/inference/forward_pass.rs`](../../packages/treetime/src/commands/timetree/inference/forward_pass.rs)
- [`packages/treetime/src/commands/timetree/utils.rs`](../../packages/treetime/src/commands/timetree/utils.rs)

Golden-master comparison of v1 timetree inference against v0 Python outputs. Runs full timetree pipeline on `flu/h3n2/20` dataset and compares node dates, branch lengths, clock rate, and coalescent contributions.

Three execution modes:

1. Poisson (parsimony branch lengths)
2. Marginal sparse (ML branch lengths, sparse sequences)
3. Marginal dense (ML branch lengths, dense sequences)

```bash
./dev/docker/run ./dev/dev E timetree_validation
```

---

## Coalescent Validation

**Test:** [`packages/treetime/examples/coalescent_validation.rs`](../../packages/treetime/examples/coalescent_validation.rs)

**Impl:** [`packages/treetime/src/commands/timetree/coalescent/coalescent.rs`](../../packages/treetime/src/commands/timetree/coalescent/coalescent.rs)

Golden-master comparison of v1 coalescent model contributions against v0 Python golden outputs. Tests piecewise constant and skyline coalescent on multiple datasets with varying Tc parameters.

```bash
./dev/docker/run ./dev/dev E coalescent_validation
```

---

## treetime-validation Crate

**Crate:** [`packages/treetime-validation/`](../../packages/treetime-validation/)

Support library for math operation validation. Provides test suites, runners, and metrics for comparing v1 numerical results against reference outputs.

### Test Suites (no `#[test]`, used by validation examples)

| Suite                            | File                                                                                                                                    | Purpose                             |
| -------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------- |
| Gaussian                         | [`gaussian.rs`](../../packages/treetime-validation/src/testing/test_suites/gaussian.rs)                                                 | Gaussian distribution operations    |
| Exponential                      | [`exponential.rs`](../../packages/treetime-validation/src/testing/test_suites/exponential.rs)                                           | Exponential distribution operations |
| Gaussian-Exponential             | [`gaussian_exponential.rs`](../../packages/treetime-validation/src/testing/test_suites/gaussian_exponential.rs)                         | Mixed distribution operations       |
| Gaussian Pairwise Multiplication | [`gaussian_pairwise_multiplication.rs`](../../packages/treetime-validation/src/testing/test_suites/gaussian_pairwise_multiplication.rs) | Pairwise product validation         |
| Gaussian Chain Multiplication    | [`gaussian_chain_multiplication.rs`](../../packages/treetime-validation/src/testing/test_suites/gaussian_chain_multiplication.rs)       | Chain product validation            |

### Metrics Tests

**Test:** [`packages/treetime-validation/src/testing/metrics/pointwise/pointwise.rs`](../../packages/treetime-validation/src/testing/metrics/pointwise/pointwise.rs) (inline)

**Impl:** [`packages/treetime-validation/src/testing/metrics/pointwise/pointwise.rs`](../../packages/treetime-validation/src/testing/metrics/pointwise/pointwise.rs)

| Test                     | Purpose                            |
| ------------------------ | ---------------------------------- |
| `test_perfect_agreement` | Zero error on identical inputs     |
| `test_constant_offset`   | Detects systematic constant offset |

**Test:** [`packages/treetime-validation/src/testing/metrics/spatial/spatial.rs`](../../packages/treetime-validation/src/testing/metrics/spatial/spatial.rs) (inline)

**Impl:** [`packages/treetime-validation/src/testing/metrics/spatial/spatial.rs`](../../packages/treetime-validation/src/testing/metrics/spatial/spatial.rs)

| Test                                    | Purpose                         |
| --------------------------------------- | ------------------------------- |
| `test_perfect_agreement`                | Zero spatial error on identical |
| `test_cumulative_error_systematic_bias` | Detects systematic spatial bias |

**Test:** [`packages/treetime-validation/src/testing/metrics/metrics.rs`](../../packages/treetime-validation/src/testing/metrics/metrics.rs) (inline)

**Impl:** [`packages/treetime-validation/src/testing/metrics/metrics.rs`](../../packages/treetime-validation/src/testing/metrics/metrics.rs)

| Test                                  | Purpose                             |
| ------------------------------------- | ----------------------------------- |
| `test_comprehensive_metrics_creation` | Default metrics config construction |
| `test_metrics_with_custom_config`     | Custom threshold metrics config     |

**Test:** [`packages/treetime-validation/src/testing/metrics/distribution/statistics.rs`](../../packages/treetime-validation/src/testing/metrics/distribution/statistics.rs) (inline)

**Impl:** [`packages/treetime-validation/src/testing/metrics/distribution/statistics.rs`](../../packages/treetime-validation/src/testing/metrics/distribution/statistics.rs)

| Test                       | Purpose                       |
| -------------------------- | ----------------------------- |
| `test_statistical_metrics` | KS statistic, mean, std tests |

**Test:** [`packages/treetime-validation/src/testing/metrics/distribution/histograms.rs`](../../packages/treetime-validation/src/testing/metrics/distribution/histograms.rs) (inline)

**Impl:** [`packages/treetime-validation/src/testing/metrics/distribution/histograms.rs`](../../packages/treetime-validation/src/testing/metrics/distribution/histograms.rs)

| Test                      | Purpose                      |
| ------------------------- | ---------------------------- |
| `test_histogram_creation` | Histogram binning and counts |
