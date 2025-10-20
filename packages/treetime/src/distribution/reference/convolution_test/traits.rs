use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;

/// Trait for function-specific convolution input generators
pub trait ConvInput: Send + Sync + Sized {
  /// Test case type for this function
  type TestCase: TestCase;

  /// Create new input with all test cases
  fn new() -> Self;

  /// Create new input with specific test cases
  fn with_test_cases(test_cases: Vec<Self::TestCase>) -> Self;

  /// Create input function f with test case parameters
  fn create_f(&self, test_case: &Self::TestCase) -> Result<GridFn, Report>;

  /// Create input function g with test case parameters
  fn create_g(&self, test_case: &Self::TestCase) -> Result<GridFn, Report>;

  /// Get evaluation domain for test case
  fn eval_domain(&self, test_case: &Self::TestCase) -> (f64, f64);

  /// Compute analytical convolution result on evaluation grid
  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<GridFn, Report>;

  /// Get all test cases
  fn test_cases(&self) -> &[Self::TestCase];

  /// Get function type name
  fn function_type(&self) -> &'static str;

  /// List all available test cases to stdout
  fn list_test_cases(&self) {
    println!("Available {} test cases:", self.function_type());
    for case in self.test_cases() {
      println!("  - {} : {}", case.name(), case.description());
    }
  }

  /// Filter test cases by comma-separated names
  fn filter_test_cases(&self, test_cases_str: &str) -> Result<Vec<Self::TestCase>, Report> {
    let all_cases = self.test_cases();
    let requested_names: Vec<&str> = test_cases_str.split(',').map(|s| s.trim()).collect();
    let filtered_cases: Vec<Self::TestCase> = all_cases
      .iter()
      .filter(|case| requested_names.contains(&case.name()))
      .cloned()
      .collect();

    if filtered_cases.is_empty() {
      eyre::bail!(
        "No matching test cases found for: {}. Available test cases: {}",
        test_cases_str,
        all_cases.iter().map(|c| c.name()).collect::<Vec<_>>().join(", ")
      );
    }

    Ok(filtered_cases)
  }
}

/// Trait for convolution algorithm implementations
pub trait ConvAlgo: Send + Sync {
  /// Compute convolution of f and g, evaluated at x_grid
  fn convolve(&self, f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report>;

  /// Get algorithm name
  fn name(&self) -> &'static str;
}
