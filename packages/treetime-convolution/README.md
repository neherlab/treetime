# treetime-convolution

Numerical convolution algorithms implementations and their testing framework for assessing correctness, accuracy and runtime performance.

## Main concepts

- **Test Suites**: Define analytical function pairs to convolve (e.g. Gaussian \* Exponential) and their analytical convolution solution (via `TestSuite` trait)

- **Algorithms**: Numerical convolution implementations (riemann, ndarray-conv, ndarray-conv-fft) (via `Algo` trait). Used to convolve function pairs defined in **Test Suites**. Numerical results are then verified against analytical solution using **Metrics**

- **Test Cases**: Individual test scenarios with specific parameters, grid configurations, and stress conditions, to evaluate numeric **Algorithm** performance (via `TestCase` trait)

- **Metrics**: aggregate (overall quality), pointwise (per-grid-point), spatial (regional patterns), distribution (statistical properties) measures for comparing numerical and analytical results

- **Runner**: Orchestrates parallel test execution, filtering, progress tracking, result collection, and reporting (console, plots, json)

## Exported Items

For use by other crates:

- `convolve()` - convolution function using ndarray-conv algorithm
- `GridFn` - function representation on a regular grid with linear interpolation

## Running Tests

### Basic usage

```bash
cargo run --example convolution_test
```

```bash
./dev/docker/run ./dev/dev E convolution_test
```

### Development mode

```bash
cargo run --example convolution_test -- --output-dir=tmp/conv-test
```

```bash
./dev/docker/run ./dev/dev E convolution_test -- --output-dir=tmp/conv-test
```

### Production mode

```bash
cargo run --release --example convolution_test -- --output-dir=tmp/conv-test
```

```bash
./dev/docker/run ./dev/dev Er convolution_test -- --output-dir=tmp/conv-test
```

### Restrict algorithms

```bash
cargo run --example convolution_test -- --algorithms=riemann,ndarray-conv
```

```bash
./dev/docker/run ./dev/dev E convolution_test -- --algorithms=riemann,ndarray-conv
```

### Restrict test suites

```bash
cargo run --example convolution_test -- --test-suites=gaussian
```

```bash
./dev/docker/run ./dev/dev E convolution_test -- --test-suites=gaussian
```

### Restrict test cases

```bash
cargo run --example convolution_test -- --test-cases=narrow,wide
```

```bash
./dev/docker/run ./dev/dev E convolution_test -- --test-cases=narrow,wide
```

### List available test cases

```bash
cargo run --example convolution_test -- --list-cases
```

```bash
./dev/docker/run ./dev/dev E convolution_test -- --list-cases
```

### Filter by slowness

```bash
cargo run --example convolution_test -- --slowness=0.3
```

```bash
./dev/docker/run ./dev/dev E convolution_test -- --slowness=0.3
```

Runs only test cases with slowness in `[0.0, 0.3]`. Higher slowness values indicate longer execution time.

### Verbose output

```bash
cargo run --example convolution_test -- --verbose
```

```bash
./dev/docker/run ./dev/dev E convolution_test -- --verbose
```

## Architecture

### Algorithms (`algos/`)

Convolution algorithm implementations:

- `riemann/` - Riemann sum approximation
- `ndarray_conv/` - Direct convolution using ndarray-conv (default, exported as `convolve()`)
- `ndarray_conv_fft/` - FFT-based convolution using ndarray-conv

All algorithms implement the `Algo` trait with `convolve()` method.

### Testing Framework (`testing/`)

- `test_suites/` - Test suite definitions (`gaussian`, `exponential`) implementing the `TestSuite` trait
- `metrics/` - Comprehensive error analysis with aggregate, pointwise, spatial, and distribution metrics
- `framework/` - Core testing infrastructure with test cases, results, and summaries
- `plots/` - Plot generation for visual analysis
- `console/` - Terminal output formatting
- `run.rs` - Main test execution logic

### Grid Functions (`grid_fn.rs`)

`GridFn` represents a function on a regular grid with:

- Linear interpolation between grid points
- Construction from arrays or domain specification
- Serialization support

## Adding Test Cases

To add test cases to an existing test suite:

1. Edit the test suite file (e.g. `src/testing/test_suites/gaussian.rs`):

   - In `GaussianTestSuite::create_test_cases()`, add new `GaussianTestCase` struct instance to the returned `vec![]`
   - Populate all struct fields.

2. Framework automatically discovers new cases through `create_test_cases()` method

## Adding Test Suites

To add a new test suite (e.g. convolving Laplacian with Gaussian):

1. Create `src/testing/test_suites/laplacian_gaussian.rs`:

   - Define test case struct with all required fields:
     ```rust
     #[derive(Debug, Clone, Serialize, Deserialize)]
     pub struct LaplacianGaussianTestCase {
       pub name: String,
       pub description: String,
       pub stress_type: String,
       pub analytical_caution: String,
       pub slowness: f64,
       pub lambda: f64,
       pub sigma: f64,
       pub mu: f64,
       pub input_grid_domain: (f64, f64),
       pub input_grid_n_points: usize,
       pub output_grid_domain: (f64, f64),
       pub output_grid_n_points: usize,
     }
     ```
   - Implement `TestCase` trait for the struct with accessor methods
   - Define suite struct with `#[derive(Default)]`:
     ```rust
     #[derive(Default)]
     pub struct LaplacianGaussianTestSuite;
     ```
   - Implement `TestSuite` trait:

     ```rust
     impl TestSuite for LaplacianGaussianTestSuite {
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

To add a convolution algorithm (e.g. Simpson's rule):

1. Create `src/algos/simpson/` directory:

   - `mod.rs` containing `pub mod simpson;`

   - `simpson.rs` containing:

     - Convolution function implementing Simpson's 1/3 rule for numerical integration:

       ```rust
       pub fn convolve_simpson(
         input_grid: &Array1<f64>,
         f_values: &Array1<f64>,
         g_values: &Array1<f64>,
         output_grid: &Array1<f64>,
       ) -> Result<Array1<f64>, Report> {
         // Read f and g over input_grid and convolve over output_grid
       }
       ```

     - Algorithm struct:

       ```rust
       pub struct SimpsonAlgo;
       ```

     - `Algo` trait implementation:

       ```rust
       impl Algo for SimpsonAlgo {
         fn name(&self) -> &'static str {
           "simpson"
         }

         fn convolve(
           &self,
           input_grid: &Array1<f64>,
           f_values: &Array1<f64>,
           g_values: &Array1<f64>,
           output_grid: &Array1<f64>,
         ) -> Result<Array1<f64>, Report> {
           convolve_simpson(input_grid, f_values, g_values, output_grid)
         }
       }
       ```

2. Edit `src/algos/mod.rs`:

   - Add `pub mod simpson;`

3. Edit `src/algos/algos.rs`:

   - Add `Simpson` variant to `ConvolutionAlgorithm` enum

   - Add match arm in `ConvolutionAlgorithm::instantiate()`:

     ```rust
     Self::Simpson => Box::new(SimpsonAlgo)
     ```

   - Add import: `use crate::algos::simpson::simpson::SimpsonAlgo;`
