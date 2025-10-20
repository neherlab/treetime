use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::grid_fn::GridFn;
use crate::make_error;
use eyre::Report;
use ndarray::Array1;
use std::collections::BTreeSet;

/// Trait for function-specific convolution input generators
pub trait ConvInput: Send + Sync {
  /// Test case type for this function
  type TestCase: TestCase;

  /// Get function type name
  fn function_type() -> &'static str;

  /// Create input function f with test case parameters
  fn create_f(test_case: &Self::TestCase) -> Result<GridFn, Report>;

  /// Create input function g with test case parameters
  fn create_g(test_case: &Self::TestCase) -> Result<GridFn, Report>;

  /// Get evaluation domain for test case
  fn eval_domain(test_case: &Self::TestCase) -> (f64, f64);

  /// Compute analytical convolution result on evaluation grid
  fn analytical_convolution(test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<GridFn, Report>;

  /// Create complete set of test cases for this function type
  fn create_test_cases() -> Vec<Self::TestCase>;

  /// List all available test cases to stdout
  fn list_test_cases() {
    println!("Available {} test cases:", Self::function_type());
    for case in Self::create_test_cases() {
      println!("  - {} : {}", case.name(), case.description());
    }
  }

  /// Create optional-filtered set of test cases
  fn filter_test_cases(filter: Option<&str>) -> Result<Vec<Self::TestCase>, Report> {
    let all_cases = Self::create_test_cases();

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
}

/// Trait for convolution algorithm implementations
pub trait ConvAlgo: Send + Sync {
  /// Compute convolution of f and g, evaluated at x_grid
  fn convolve(&self, f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report>;

  /// Get algorithm name
  fn name(&self) -> &'static str;
}
