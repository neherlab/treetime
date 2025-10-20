use crate::distribution::reference::convolution_test::algorithms::ConvolutionAlgorithm;
use crate::distribution::reference::convolution_test::framework::results::TestResult;
use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::convolution_test::metrics::metrics::ConvolutionMetrics;
use crate::distribution::reference::grid_fn::GridFn;
use crate::make_error;
use eyre::Report;
use ndarray::Array1;
use std::collections::BTreeSet;
use std::time::Instant;

/// Trait for function-specific convolution input generators and test runners
pub trait ConvInput: Send + Sync {
  /// Test case type for this function
  type TestCase: TestCase;

  /// Get function type name
  fn function_type(&self) -> &'static str;

  /// Create input function f with test case parameters
  fn create_f(&self, test_case: &Self::TestCase) -> Result<GridFn, Report>;

  /// Create input function g with test case parameters
  fn create_g(&self, test_case: &Self::TestCase) -> Result<GridFn, Report>;

  /// Get evaluation domain for test case
  fn eval_domain(&self, test_case: &Self::TestCase) -> (f64, f64);

  /// Compute analytical convolution result on evaluation grid
  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<GridFn, Report>;

  /// Create complete set of test cases for this function type
  fn create_test_cases(&self) -> Vec<Self::TestCase>;

  /// List all available test cases to stdout
  fn list_test_cases(&self) {
    println!("Available {} test cases:", self.function_type());
    for case in self.create_test_cases() {
      println!("  - {} : {}", case.name(), case.description());
    }
  }

  /// Create optional-filtered set of test cases
  fn filter_test_cases(&self, filter: Option<&str>) -> Result<Vec<Self::TestCase>, Report> {
    let all_cases = self.create_test_cases();

    let filter = filter.and_then(|value| {
      let trimmed = value.trim();
      if trimmed.is_empty() { None } else { Some(trimmed) }
    });

    match filter {
      None | Some("all") => Ok(all_cases),
      Some(filter_str) => {
        let requested_names: BTreeSet<&str> = filter_str.split(',').map(|name| name.trim()).collect();
        let filtered_cases: Vec<Self::TestCase> = all_cases
          .iter()
          .filter(|case| requested_names.contains(case.name()))
          .cloned()
          .collect();

        if filtered_cases.is_empty() {
          let available_names = all_cases.iter().map(|case| case.name()).collect::<Vec<_>>().join(", ");
          return make_error!(
            "No matching test cases found for: {filter_str}. Available test cases: {available_names}"
          );
        }

        Ok(filtered_cases)
      },
    }
  }

  /// Run a test for the given test case and algorithm
  fn run_test(
    &self,
    test_case: &Self::TestCase,
    algorithm: ConvolutionAlgorithm,
  ) -> Result<TestResult<Self::TestCase>, Report> {
    let start_time = Instant::now();

    let f = self.create_f(test_case)?;
    let g = self.create_g(test_case)?;

    let (eval_min, eval_max) = self.eval_domain(test_case);
    let n_eval_points = ((eval_max - eval_min) / test_case.dx() + 1.0).round() as usize;
    let eval_grid = Array1::from_iter((0..n_eval_points).map(|i| eval_min + i as f64 * test_case.dx()));

    let algo: Box<dyn ConvAlgo> = match algorithm {
      ConvolutionAlgorithm::Riemann => {
        Box::new(crate::distribution::reference::convolution_test::algo_impls::RiemannAlgo)
      },
      ConvolutionAlgorithm::Ndarray => {
        Box::new(crate::distribution::reference::convolution_test::algo_impls::NdarrayAlgo)
      },
    };

    let actual_result = algo.convolve(&f, &g, &eval_grid)?;
    let expected_result = self.analytical_convolution(test_case, &eval_grid)?;

    let execution_time = start_time.elapsed().as_secs_f64() * 1000.0;

    let metrics = ConvolutionMetrics::new(
      actual_result.x(),
      actual_result.y(),
      expected_result.y(),
      execution_time,
    )?;

    let evaluation_grid = actual_result.x().to_owned();
    let actual_values = actual_result.y().to_owned();
    let expected_values = expected_result.y().to_owned();
    let f_x_values = f.x().to_owned();
    let f_y_values = f.y().to_owned();
    let g_x_values = g.x().to_owned();
    let g_y_values = g.y().to_owned();

    Ok(TestResult {
      algorithm,
      test_case: test_case.clone(),
      execution_time_ms: execution_time,
      f_x_values,
      f_y_values,
      g_x_values,
      g_y_values,
      evaluation_grid,
      actual_values,
      expected_values,
      metrics,
    })
  }
}

/// Trait for convolution algorithm implementations
pub trait ConvAlgo: Send + Sync {
  /// Compute convolution of f and g, evaluated at x_grid
  fn convolve(&self, f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report>;

  /// Get algorithm name
  fn name(&self) -> &'static str;
}
