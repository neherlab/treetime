# Validation & Smoke Tests

[Back to index](_index.md)

## Summary

| Category                  | Tests   | Type          |
| ------------------------- | ------- | ------------- |
| Smoke tests (CLI)         | 175     | Integration   |
| Timetree validation       | 1       | Golden-master |
| Coalescent validation     | 1       | Golden-master |
| treetime-validation crate | 8       | Unit          |
| **Total**                 | **185** |               |

---

## Smoke Tests

**Script:** [`dev/run-smoke-tests`](../../dev/run-smoke-tests)

End-to-end CLI validation running the treetime binary across datasets and flag combinations. Executed in parallel via Python multiprocessing.

### Breakdown

| Command   | Datasets | Modes                                      | Runs    |
| --------- | -------- | ------------------------------------------ | ------- |
| ancestral | 17       | parsimony, marginal-sparse, marginal-dense | 51      |
| clock     | 17       | default                                    | 17      |
| optimize  | 17       | sparse (dense disabled)                    | 17      |
| prune     | 17       | default                                    | 17      |
| timetree  | 67       | varied flag combinations                   | 67      |
| help      | 6        | --help per subcommand                      | 6       |
| **Total** |          |                                            | **175** |

### Datasets Used

`dengue/{20,100}`, `ebola/{20,100,362}`, `flu/h3n2/{20,200}`, `lassa/L/{20,50,200}`, `rsv/a/{20,100}`, `tb/{20,100,149}`, `zika/{20,86}`

Additional timetree-only: `mpox/clade-ii/20`

### Timetree Flag Coverage (67 runs)

| Category           | Flags tested                                                    | Runs |
| ------------------ | --------------------------------------------------------------- | ---- |
| Basic              | defaults                                                        | 5    |
| Clock rate         | `--clock-rate=0.003`                                            | 2    |
| Branch length mode | `--branch-length-mode=input`                                    | 2    |
| Time marginal      | `--time-marginal={always,only-final}`, `--confidence`           | 6    |
| Polytomies         | `--keep-polytomies`, `--resolve-polytomies`                     | 3    |
| Coalescent         | `--coalescent=10.0`, `--coalescent-opt`, `--coalescent-skyline` | 11   |
| Reroot             | `--reroot={min-dev,oldest,least-squares}`, `--keep-root`        | 7    |
| Clock filter       | `--clock-filter={0,5.0}`                                        | 2    |
| Covariation        | `--covariation`                                                 | 2    |
| Ancestral method   | `--method-anc=parsimony`                                        | 1    |
| GTR model          | `--gtr=jc69`                                                    | 1    |
| Max iterations     | `--max-iter={0,1,4}`                                            | 4    |
| Tip slack          | `--tip-slack=0.1`                                               | 1    |
| Combinations       | multi-flag combinations                                         | 9    |
| Known failures     | `--vary-rate`, `--relax`, `--dense=true`, crash-prone datasets  | 12   |

---

## Timetree Validation

**File:** [`examples/timetree_validation.rs`](../../packages/treetime/examples/timetree_validation.rs)

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

**File:** [`examples/coalescent_validation.rs`](../../packages/treetime/examples/coalescent_validation.rs)

Golden-master comparison of v1 coalescent model contributions against v0 Python golden outputs. Tests piecewise constant and skyline coalescent on multiple datasets with varying Tc parameters.

```bash
./dev/docker/run ./dev/dev E coalescent_validation
```

---

## treetime-validation Crate

**Crate:** [`packages/treetime-validation/`](../../packages/treetime-validation/)

Support library for math operation validation. Provides test suites, runners, and metrics for comparing v1 numerical results against reference outputs.

### Test Suites (no `#[test]`, used by validation examples)

| Suite                            | File                                              | Purpose                             |
| -------------------------------- | ------------------------------------------------- | ----------------------------------- |
| Gaussian                         | `test_suites/gaussian.rs`                         | Gaussian distribution operations    |
| Exponential                      | `test_suites/exponential.rs`                      | Exponential distribution operations |
| Gaussian-Exponential             | `test_suites/gaussian_exponential.rs`             | Mixed distribution operations       |
| Gaussian Pairwise Multiplication | `test_suites/gaussian_pairwise_multiplication.rs` | Pairwise product validation         |
| Gaussian Chain Multiplication    | `test_suites/gaussian_chain_multiplication.rs`    | Chain product validation            |

### Metrics Tests (8 `#[test]` functions)

**File:** `testing/metrics/pointwise/pointwise.rs` (inline)

| Test                     | Purpose                            |
| ------------------------ | ---------------------------------- |
| `test_perfect_agreement` | Zero error on identical inputs     |
| `test_constant_offset`   | Detects systematic constant offset |

**File:** `testing/metrics/spatial/spatial.rs` (inline)

| Test                                    | Purpose                         |
| --------------------------------------- | ------------------------------- |
| `test_perfect_agreement`                | Zero spatial error on identical |
| `test_cumulative_error_systematic_bias` | Detects systematic spatial bias |

**File:** `testing/metrics/metrics.rs` (inline)

| Test                                  | Purpose                             |
| ------------------------------------- | ----------------------------------- |
| `test_comprehensive_metrics_creation` | Default metrics config construction |
| `test_metrics_with_custom_config`     | Custom threshold metrics config     |

**File:** `testing/metrics/distribution/statistics.rs` (inline)

| Test                       | Purpose                       |
| -------------------------- | ----------------------------- |
| `test_statistical_metrics` | KS statistic, mean, std tests |

**File:** `testing/metrics/distribution/histograms.rs` (inline)

| Test                      | Purpose                      |
| ------------------------- | ---------------------------- |
| `test_histogram_creation` | Histogram binning and counts |
