#![allow(clippy::many_single_char_names)]
use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::convolution_test::traits::ConvInput;
use crate::distribution::reference::grid_fn::GridFn;
use crate::make_error;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub struct ExponentialConvInput {
  test_cases: Vec<ExponentialTestCase>,
}

impl ExponentialConvInput {
  pub fn new(test_case_filter: &str) -> Result<Self, Report> {
    let all_cases = Self::create_all_test_cases();
    let test_cases = Self::filter_test_cases_impl(&all_cases, Some(test_case_filter))?;
    Ok(Self { test_cases })
  }

  fn create_all_test_cases() -> Vec<ExponentialTestCase> {
    vec![
      ExponentialTestCase {
        name: "moderate_coarse_grid".to_owned(),
        description: "Moderate rates on a coarse grid. Discretization error on coarse grids; step size sensitivity."
          .to_owned(),
        stress_type: "aliasing/accuracy vs Δx, Δx scaling factor".to_owned(),
        analytical_caution: "none".to_owned(),
        a: 1.0,
        b: 0.8,
        f_domain: (0.0, 20.0),
        g_domain: (0.0, 20.0),
        eval_domain: (0.0, 40.0),
        dx: 0.1,
      },
      ExponentialTestCase {
        name: "tight_truncation".to_owned(),
        description:
          "Tight truncation to expose intentional errors. Sensitivity to insufficient support; boundary effects."
            .to_owned(),
        stress_type: "wrap-around artifacts if padding insufficient".to_owned(),
        analytical_caution: "none".to_owned(),
        a: 1.0,
        b: 2.0,
        f_domain: (0.0, 3.0),
        g_domain: (0.0, 2.0),
        eval_domain: (0.0, 5.0),
        dx: 0.01,
      },
    ]
  }

  fn filter_test_cases_impl(all_cases: &[ExponentialTestCase], filter: Option<&str>) -> Result<Vec<ExponentialTestCase>, Report> {
    let filter = filter.and_then(|value| {
      let trimmed = value.trim();
      if trimmed.is_empty() { None } else { Some(trimmed) }
    });

    match filter {
      None | Some("all") => Ok(all_cases.to_vec()),
      Some(filter_str) => {
        let requested_names: BTreeSet<&str> = filter_str.split(',').map(|name| name.trim()).collect();
        let filtered_cases: Vec<ExponentialTestCase> = all_cases
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

impl ConvInput for ExponentialConvInput {
  type TestCase = ExponentialTestCase;

  fn function_type(&self) -> &'static str {
    "exponential"
  }

  fn test_cases(&self) -> &[Self::TestCase] {
    &self.test_cases
  }

  fn create_f(&self, test_case: &Self::TestCase) -> Result<GridFn, Report> {
    let ExponentialTestCase { a, f_domain, dx, .. } = test_case;
    // Generate exponential function f(x) = a*exp(-ax) for x ≥ 0, 0 otherwise
    GridFn::from_grid(*f_domain, *dx, |x| if x >= 0.0 { a * (-a * x).exp() } else { 0.0 })
  }

  fn create_g(&self, test_case: &Self::TestCase) -> Result<GridFn, Report> {
    let ExponentialTestCase { b, g_domain, dx, .. } = test_case;
    // Generate exponential function g(x) = b*exp(-bx) for x ≥ 0, 0 otherwise
    GridFn::from_grid(*g_domain, *dx, |x| if x >= 0.0 { b * (-b * x).exp() } else { 0.0 })
  }

  fn eval_domain(&self, test_case: &Self::TestCase) -> (f64, f64) {
    test_case.eval_domain
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, _eval_grid: &Array1<f64>) -> Result<GridFn, Report> {
    let ExponentialTestCase {
      a, b, eval_domain, dx, ..
    } = test_case;
    // Analytical convolution of two exponential functions
    //
    // f(x) = a*exp(-ax) for x ≥ 0, 0 otherwise
    // g(x) = b*exp(-bx) for x ≥ 0, 0 otherwise
    //
    // Result: (ab)/(a-b) * (1-exp(-(a-b)x)) * exp(-bx) for a ≠ b
    // Special case: when a = b, result is ab*x*exp(-ax)
    GridFn::from_grid(*eval_domain, *dx, |x| {
      if x < 0.0 {
        0.0
      } else if (a - b).abs() < 1e-15 {
        a * b * x * (-a * x).exp()
      } else {
        (a * b) / (a - b) * (1.0 - (-(a - b) * x).exp()) * (-b * x).exp()
      }
    })
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    Self::create_all_test_cases()
  }

  fn filter_test_cases(&self, filter: Option<&str>) -> Result<Vec<Self::TestCase>, Report> {
    let all_cases = self.create_test_cases();
    Self::filter_test_cases_impl(&all_cases, filter)
  }
}

/// Test case parameters for Exponential convolution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExponentialTestCase {
  pub name: String,
  pub description: String,
  pub stress_type: String,
  pub analytical_caution: String,
  pub a: f64,                  // Decay rate for f(x) = exp(-ax) for x >= 0
  pub b: f64,                  // Decay rate for g(x) = exp(-bx) for x >= 0
  pub f_domain: (f64, f64),    // Should start at 0 for causal support
  pub g_domain: (f64, f64),    // Should start at 0 for causal support
  pub eval_domain: (f64, f64), // Should start at 0 for causal support
  pub dx: f64,
}

impl TestCase for ExponentialTestCase {
  fn name(&self) -> &str {
    &self.name
  }

  fn description(&self) -> &str {
    &self.description
  }

  fn stress_type(&self) -> &str {
    &self.stress_type
  }

  fn analytical_caution(&self) -> &str {
    &self.analytical_caution
  }

  fn dx(&self) -> f64 {
    self.dx
  }
}
