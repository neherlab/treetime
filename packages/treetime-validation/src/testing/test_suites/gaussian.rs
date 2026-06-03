use crate::testing::framework::test_case::{TestCase, TestCaseBase};
use crate::testing::test_suites::test_suites::ConvolutionTestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use treetime_analytical::validation::cases::gaussian_convolution::{
  GAUSSIAN_CONVOLUTION_CASES, GaussianConvolutionTestCase,
};
use treetime_analytical::{gaussian_convolution_pdf_grid, gaussian_pdf_grid};

#[derive(Default)]
pub struct GaussianTestSuite;

impl ConvolutionTestSuite for GaussianTestSuite {
  type TestCase = GaussianTestCase;

  fn test_suite_name(&self) -> &'static str {
    "conv-gaussian-gaussian"
  }

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(gaussian_pdf_grid(0.0, test_case.sigma_f, grid))
  }

  fn create_g(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(gaussian_pdf_grid(test_case.mu, test_case.sigma_g, grid))
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(gaussian_convolution_pdf_grid(
      test_case.sigma_f,
      test_case.sigma_g,
      test_case.mu,
      eval_grid,
    ))
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    GAUSSIAN_CONVOLUTION_CASES.iter().map(GaussianTestCase::from).collect()
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianTestCase {
  #[serde(flatten)]
  pub base: TestCaseBase,
  pub sigma_f: f64,
  pub sigma_g: f64,
  pub mu: f64,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
}

impl From<&GaussianConvolutionTestCase> for GaussianTestCase {
  fn from(case: &GaussianConvolutionTestCase) -> Self {
    Self {
      base: TestCaseBase {
        name: case.name.to_owned(),
        description: case.description.to_owned(),
        stress_type: case.stress_type.to_owned(),
        analytical_caution: case.analytical_caution.to_owned(),
        slowness: case.slowness,
      },
      sigma_f: case.sigma_f,
      sigma_g: case.sigma_g,
      mu: case.mu,
      input_grid_domain: case.input_grid_domain,
      input_grid_n_points: case.input_grid_n_points,
    }
  }
}

impl TestCase for GaussianTestCase {
  fn base(&self) -> &TestCaseBase {
    &self.base
  }

  fn input_grid_domain(&self) -> (f64, f64) {
    self.input_grid_domain
  }

  fn input_grid_n_points(&self) -> usize {
    self.input_grid_n_points
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  fn sample() -> GaussianTestCase {
    GaussianTestCase {
      base: TestCaseBase {
        name: "case-1".to_owned(),
        description: "desc".to_owned(),
        stress_type: "stress".to_owned(),
        analytical_caution: "caution".to_owned(),
        slowness: 1.5,
      },
      sigma_f: 1.0,
      sigma_g: 2.0,
      mu: 0.5,
      input_grid_domain: (-5.0, 5.0),
      input_grid_n_points: 100,
    }
  }

  // The trait accessors are defaulted via `self.base()`; verify they forward
  // the embedded `TestCaseBase` values rather than returning anything else.
  #[test]
  fn test_gaussian_test_case_accessors_forward_to_base() {
    let case = sample();
    assert_eq!("case-1", case.name());
    assert_eq!("desc", case.description());
    assert_eq!("stress", case.stress_type());
    assert_eq!("caution", case.analytical_caution());
    #[allow(clippy::float_cmp)] // 1.5 is exactly representable
    {
      assert_eq!(1.5, case.slowness());
    }
  }

  // `#[serde(flatten)]` keeps the metadata fields at the top level. A
  // regression that drops the attribute would nest them under a `base` key and
  // change the JSON written to result files and the console.
  #[test]
  fn test_gaussian_test_case_serializes_flat() -> Result<(), Report> {
    let actual = json_write_str(&sample(), JsonPretty(true))?;
    let expected = r#"{
  "name": "case-1",
  "description": "desc",
  "stress_type": "stress",
  "analytical_caution": "caution",
  "slowness": 1.5,
  "sigma_f": 1.0,
  "sigma_g": 2.0,
  "mu": 0.5,
  "input_grid_domain": [
    -5.0,
    5.0
  ],
  "input_grid_n_points": 100
}"#;
    assert_eq!(expected, actual);
    Ok(())
  }
}
