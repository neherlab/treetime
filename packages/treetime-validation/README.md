# treetime-validation

Validation framework for numerical convolution and multiplication algorithms, assessing correctness, accuracy and runtime performance.

## Main concepts

- **Test Suites**: Define analytical function pairs to convolve or multiply (e.g. Gaussian \* Exponential) and their analytical solution (via `ConvolutionTestSuite`, `MultiplicationTestSuite`, `ChainMultiplicationTestSuite` traits)

- **Algorithms**: Numerical convolution and multiplication implementations (riemann, ndarray-conv, ndarray-conv-fft, naive-multiplication, log-scale-multiplication, aggressive-multiplication). Numerical results are verified against analytical solutions using **Metrics**

- **Test Cases**: Individual test scenarios with specific parameters, grid configurations, and stress conditions, to evaluate numeric **Algorithm** performance (via `TestCase` trait)

- **Metrics**: aggregate (overall quality), pointwise (per-grid-point), spatial (regional patterns), distribution (statistical properties) measures for comparing numerical and analytical results

- **Runner**: Orchestrates parallel test execution, filtering, progress tracking, result collection, and reporting (console, plots, json)

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

### Restrict algorithms

Run all algorithms (default):

```bash
cargo run --example validation_test -- --algorithms=all
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --algorithms=all
```

Run specific algorithms (comma-separated):

```bash
cargo run --example validation_test -- --algorithms=riemann,ndarray-conv
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --algorithms=riemann,ndarray-conv
```

### Restrict multiplication algorithms

Run specific multiplication algorithms (comma-separated):

```bash
cargo run --example validation_test -- --mult-algorithms=naive-multiplication,log-scale-multiplication
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --mult-algorithms=naive-multiplication,log-scale-multiplication
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
cargo run --example validation_test -- --test-suites=gaussian,exponential
```

```bash
./dev/docker/run ./dev/dev E validation_test -- --test-suites=gaussian,exponential
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
- `MultiplicationAlgorithm` - NaiveMultiplication, LogScaleMultiplication, AggressiveMultiplication

Actual implementations live in the `treetime-ops` crate.

### Testing Framework (`testing/`)

- `test_suites/` - Test suite definitions (`gaussian`, `exponential`, `gaussian_multiplication`, etc.) implementing suite traits
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
         "laplacian-gaussian"
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
   - Add `LaplacianGaussian` variant to `TestSuiteName` enum

   - Add match arm in `TestSuiteName::run_tests()`:

     ```rust
     Self::LaplacianGaussian => run_convolution_tests_impl::<LaplacianGaussianTestSuite>(args)
     ```

   - Add import: `use crate::testing::test_suites::laplacian_gaussian::LaplacianGaussianTestSuite;`

## Adding Algorithms

Convolution and multiplication algorithms are implemented in the `treetime-ops` crate. To add a new algorithm:

1. Implement the algorithm in `treetime-ops`

2. Edit `src/algorithms.rs` in this crate:
   - Add variant to `ConvolutionAlgorithm` or `MultiplicationAlgorithm` enum

   - Add match arm in `instantiate()` method to create the algorithm instance
