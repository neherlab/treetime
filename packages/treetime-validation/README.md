# treetime-validation

Validation framework for numerical convolution and multiplication algorithms, assessing correctness, accuracy and runtime performance.

## Goals

Validate numerical convolution and multiplication algorithms used in treetime's timetree inference. These operations combine probability distributions during message passing on phylogenetic trees:

- Convolution: Propagates time distributions along branches (child time + branch length)
- Multiplication: Combines constraints at nodes (intersection of probability beliefs)

The framework measures how closely discrete numerical results match analytical (continuous) ground truth across varied function pairs, grid configurations, and stress conditions.

### What the framework validates

1. Correctness: Numerical algorithms produce mathematically correct results
2. Accuracy: Discretization error is within acceptable bounds for production use
3. Performance: Execution time comparison across algorithm variants
4. Robustness: Behavior under stress conditions (coarse grids, tight truncation, extreme parameters)

### Interpreting results

Convolution: All convolution algorithms (Riemann, NdarrayConv, NdarrayConvFft) compute the same discrete convolution sum. They produce identical accuracy metrics because they implement the same mathematical operation - differences appear only in execution time and floating-point rounding (typically 1e-15 level). Accuracy variation across test cases reflects discretization error (grid resolution, domain coverage), not algorithm quality.

Multiplication: Pairwise multiplication algorithms use identical code paths (simple element-wise product). Differences emerge only in chain multiplication (`multiply_many`) where normalization strategies differ:
- Pointwise: No intermediate normalization - can underflow to zero
- LogScale: Normalizes only when max falls below 1e-100
- Aggressive: Normalizes after every multiplication step

## Challenges

### Multiplication stress conditions

Underflow:
- Many small probability factors (deep trees, many constraints)
- Tail regions with values approaching machine epsilon
- Coalescent contributions with very small survival probabilities

Dynamic range:
- Peak values ~1.0 alongside tail values ~1e-300
- Precision loss in small values when operating with large ones
- Log-scale representation helps but requires conversion overhead

Accumulation:
- Rounding errors compound over chain length
- Normalization frequency affects accumulated error
- Order of operations can affect final precision

### Convolution stress conditions

Discretization:
- Coarse grids miss fine structure
- Grid spacing must resolve narrowest distribution
- Step size affects integral approximation accuracy

Domain truncation:
- Cutting off distribution tails loses probability mass
- Tight domains cause boundary artifacts
- Wide distributions need larger domains

Scale mismatch:
- Convolving narrow (delta-like) with wide distributions
- Very different sigma values stress grid requirements
- Result domain is sum of input domains

Translation:
- Large shifts test domain coverage
- Sub-grid shifts affect interpolation accuracy
- Negative shifts (backward pass) vs positive shifts (forward pass)

Tail precision:
- Small tail contributions must accumulate correctly
- Floating-point addition of many tiny values
- Tail region errors often larger than peak region

Grid alignment:
- Input grids with different spacing require resampling
- Resampling introduces interpolation error
- FFT requires uniform grid spacing

## Main concepts

- Test Suites: Define analytical function pairs to convolve or multiply (e.g. Gaussian \* Exponential) and their analytical solution (via `ConvolutionTestSuite`, `MultiplicationTestSuite`, `ChainMultiplicationTestSuite` traits)
- Algorithms: Numerical convolution and multiplication implementations (riemann, ndarray-conv, ndarray-conv-fft, pointwise, log-scale, aggressive). Numerical results are verified against analytical solutions using Metrics
- Test Cases: Individual test scenarios with specific parameters, grid configurations, and stress conditions, to evaluate numeric Algorithm performance (via `TestCase` trait)
- Metrics: aggregate (overall quality), pointwise (per-grid-point), spatial (regional patterns), distribution (statistical properties) measures for comparing numerical and analytical results
- Runner: Orchestrates parallel test execution, filtering, progress tracking, result collection, and reporting (console, plots, json)

## Running Tests

### Basic usage

```bash
cargo run --example validation_test
```

```bash
./dev/docker/run ./dev/dev E validation_test
```

### Development mode

```bash
cargo run --example validation_test -- --output-dir=tmp/validation-test
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --output-dir=tmp/validation-test
```

### Production mode

```bash
cargo run --release --example validation_test -- --output-dir=tmp/validation-test
```

```bash
./dev/docker/run ./dev/dev Er validation_test -- --output-dir=tmp/validation-test
```

### Restrict convolution algorithms

Run all algorithms (default):

```bash
cargo run --example validation_test -- --conv-algorithms=all
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --conv-algorithms=all
```

Run specific algorithms (comma-separated):

```bash
cargo run --example validation_test -- --conv-algorithms=riemann,ndarray-conv
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --conv-algorithms=riemann,ndarray-conv
```

### Restrict multiplication algorithms

Run specific multiplication algorithms (comma-separated):

```bash
cargo run --example validation_test -- --mult-algorithms=pointwise,log-scale
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --mult-algorithms=pointwise,log-scale
```

### Restrict test suites

Run all test suites (default):

```bash
cargo run --example validation_test -- --test-suites=all
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --test-suites=all
```

Run specific test suites (comma-separated):

```bash
cargo run --example validation_test -- --test-suites=conv-gaussian-gaussian,conv-exponential-exponential
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --test-suites=conv-gaussian-gaussian,conv-exponential-exponential
```

### Restrict test cases

Run specific test cases (comma-separated):

```bash
cargo run --example validation_test -- --test-cases=python_notebook_case_1,coarse_grid
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --test-cases=python_notebook_case_1,coarse_grid
```

### List available test cases

```bash
cargo run --example validation_test -- --list-cases
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --list-cases
```

### Filter by slowness

```bash
cargo run --example validation_test -- --slowness=0.3
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --slowness=0.3
```

Runs only test cases with slowness in `[0.0, 0.3]`. Higher slowness values indicate longer execution time.

### Verbose output

```bash
cargo run --example validation_test -- --verbose
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --verbose
```

## Architecture

### Algorithms (`algorithms.rs`)

Algorithm enums for CLI selection:

- `ConvolutionAlgorithm` - Riemann, NdarrayConv, NdarrayConvFft
- `MultiplicationAlgorithm` - Pointwise, LogScale, Aggressive

Actual implementations live in the `treetime-ops` crate.

### Testing Framework (`testing/`)

- `test_suites/` - Test suite definitions (`conv-gaussian-gaussian`, `conv-exponential-exponential`, `mult-gaussian-pairwise`, etc.) implementing suite traits
- `metrics/` - Comprehensive error analysis with aggregate, pointwise, spatial, and distribution metrics
- `framework/` - Core testing infrastructure with test cases, results, and summaries
- `plots/` - Plot generation for visual analysis
- `console/` - Terminal output formatting
- `runners/` - Test execution for convolution, multiplication, and chain multiplication
- `run.rs` - Main test execution logic

## Adding Test Cases

To add test cases to an existing test suite:

1. Edit the test suite file (e.g. `src/testing/test_suites/gaussian.rs`):
   - In `GaussianTestSuite::create_test_cases()`, add new `GaussianTestCase` struct instance to the returned `vec![]`
   - Populate all struct fields.

2. Framework automatically discovers new cases through `create_test_cases()` method

## Adding Test Suites

To add a new test suite (e.g. convolving Laplacian with Gaussian):

1. Create `src/testing/test_suites/laplacian_gaussian.rs`:
   - Define test case struct with required `TestCase` trait fields plus suite-specific parameters:
     ```rust
     #[derive(Debug, Clone, Serialize, Deserialize)]
     pub struct LaplacianGaussianTestCase {
       pub name: String,
       pub description: String,
       pub stress_type: String,
       pub analytical_caution: String,
       pub slowness: f64,
       pub input_grid_domain: (f64, f64),
       pub input_grid_n_points: usize,
       // Suite-specific parameters
       pub lambda: f64,
       pub sigma: f64,
       pub mu: f64,
     }
     ```
   - Implement `TestCase` trait for the struct with accessor methods
   - Define suite struct with `#[derive(Default)]`:
     ```rust
     #[derive(Default)]
     pub struct LaplacianGaussianTestSuite;
     ```
   - Implement `ConvolutionTestSuite` trait:

     ```rust
     impl ConvolutionTestSuite for LaplacianGaussianTestSuite {
       type TestCase = LaplacianGaussianTestCase;

       fn test_suite_name(&self) -> &'static str {
         "conv-laplacian-gaussian"
       }

       fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
         let LaplacianGaussianTestCase { lambda, .. } = test_case;
         // Compute Laplacian on grid using lambda parameter
       }

       fn create_g(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
         let LaplacianGaussianTestCase { sigma, mu, .. } = test_case;
         // Compute Gaussian on grid using sigma and mu parameters
       }

       fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
         let LaplacianGaussianTestCase { lambda, sigma, mu, .. } = test_case;
         // Compute analytical Laplacian * Gaussian convolution on eval_grid using lambda, sigma, mu
       }

       fn create_test_cases(&self) -> Vec<Self::TestCase> {
         vec![
           LaplacianGaussianTestCase {
             name: "baseline".to_owned(),
             // ... populate all fields
           },
           // Add more test cases
         ]
       }
     }
     ```

2. Edit `src/testing/test_suites/mod.rs`:
   - Add `pub mod laplacian_gaussian;`

3. Edit `src/testing/test_suites/test_suites.rs`:
   - Add `ConvLaplacianGaussian` variant to `TestSuiteName` enum

   - Add match arm in `TestSuiteName::run_tests()`:

     ```rust
     Self::ConvLaplacianGaussian => run_convolution_tests_impl::<LaplacianGaussianTestSuite>(args)
     ```

   - Add import: `use crate::testing::test_suites::laplacian_gaussian::LaplacianGaussianTestSuite;`

## Adding Algorithms

Convolution and multiplication algorithms are implemented in the `treetime-ops` crate. To add a new algorithm:

1. Implement the algorithm in `treetime-ops`

2. Edit `src/algorithms.rs` in this crate:
   - Add variant to `ConvolutionAlgorithm` or `MultiplicationAlgorithm` enum

   - Add match arm in `instantiate()` method to create the algorithm instance
