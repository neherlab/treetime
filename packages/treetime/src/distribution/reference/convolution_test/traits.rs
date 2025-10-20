use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;

/// Trait for function-specific convolution input generators
pub trait ConvInput: Send + Sync + Sized {
  /// Test case type for this function
  type TestCase: TestCase;

  /// Get function type name
  fn function_type(&self) -> &'static str;

  /// Create new input with specific test cases
  fn new_with_cases(test_cases: Vec<Self::TestCase>) -> Self;

  /// Create new input with all test cases
  fn new_all() -> Self {
    Self::new_with_cases(Self::create_test_cases())
  }

  /// Create new input with optional test case filter
  fn new(filter: Option<&str>) -> Result<Self, Report> {
    match filter {
      None | Some("all") => Ok(Self::new_all()),
      Some(filter_str) => {
        let all_cases = Self::create_test_cases();
        let requested_names: Vec<&str> = filter_str.split(',').map(|s| s.trim()).collect();
        let filtered_cases: Vec<Self::TestCase> = all_cases
          .iter()
          .filter(|case| requested_names.contains(&case.name()))
          .cloned()
          .collect();

        if filtered_cases.is_empty() {
          eyre::bail!(
            "No matching test cases found for: {}. Available test cases: {}",
            filter_str,
            all_cases.iter().map(|c| c.name()).collect::<Vec<_>>().join(", ")
          );
        }

        Ok(Self::new_with_cases(filtered_cases))
      },
    }
  }

  /// Create input function f with test case parameters
  fn create_f(&self, test_case: &Self::TestCase) -> Result<GridFn, Report>;

  /// Create input function g with test case parameters
  fn create_g(&self, test_case: &Self::TestCase) -> Result<GridFn, Report>;

  /// Get evaluation domain for test case
  fn eval_domain(&self, test_case: &Self::TestCase) -> (f64, f64);

  /// Compute analytical convolution result on evaluation grid
  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<GridFn, Report>;

  /// List all available test cases to stdout
  fn list_test_cases(&self) {
    println!("Available {} test cases:", self.function_type());
    for case in self.test_cases() {
      println!("  - {} : {}", case.name(), case.description());
    }
  }

  /// Get all test cases
  fn test_cases(&self) -> &[Self::TestCase];

  /// Create complete set of test cases for this function type
  fn create_test_cases() -> Vec<Self::TestCase>;
}

/// Trait for convolution algorithm implementations
pub trait ConvAlgo: Send + Sync {
  /// Compute convolution of f and g, evaluated at x_grid
  fn convolve(&self, f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report>;

  /// Get algorithm name
  fn name(&self) -> &'static str;
}
