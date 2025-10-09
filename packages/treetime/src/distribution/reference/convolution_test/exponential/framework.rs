use crate::distribution::reference::convolution::{ndarray_convolve, riemann_convolve};
use crate::distribution::reference::convolution_test::algorithms::ConvolutionAlgorithm;
use crate::distribution::reference::convolution_test::framework::{ConvolutionTestRunner, TestResult};
use crate::distribution::reference::convolution_test::output::ToFlatResult;
use crate::distribution::reference::domain_agreement_metrics::DomainAgreementMetrics;
use crate::distribution::reference::exponential::{exponential_convolution, exponential_f, exponential_g};
use eyre::Report;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use std::time::Instant;

use super::test_cases::{ExponentialTestCase, create_exponential_test_cases};

/// Exponential-specific test framework implementation
#[derive(Clone, Debug)]
pub struct ExponentialTestRunner {
  test_cases: Vec<ExponentialTestCase>,
}

impl ExponentialTestRunner {
  /// Create new Exponential test runner with default test cases
  pub fn new() -> Self {
    Self {
      test_cases: create_exponential_test_cases(),
    }
  }

  /// Create with custom test cases
  pub fn with_test_cases(test_cases: Vec<ExponentialTestCase>) -> Self {
    Self { test_cases }
  }
}

impl ConvolutionTestRunner<ExponentialTestCase> for ExponentialTestRunner {
  fn run_test(
    &self,
    test_case: &ExponentialTestCase,
    algorithm: ConvolutionAlgorithm,
  ) -> Result<TestResult<ExponentialTestCase>, Report> {
    let start_time = Instant::now();

    // Create input functions
    let f = exponential_f(test_case.a, test_case.f_domain, test_case.dx)?;
    let g = exponential_g(test_case.b, test_case.g_domain, test_case.dx)?;

    // Create evaluation grid
    let (eval_min, eval_max) = test_case.eval_domain;
    let n_eval_points = ((eval_max - eval_min) / test_case.dx + 1.0).round() as usize;
    let eval_grid = Array1::from_iter((0..n_eval_points).map(|i| eval_min + i as f64 * test_case.dx));

    // Run convolution
    let actual_result = match algorithm {
      ConvolutionAlgorithm::Riemann => riemann_convolve(&f, &g, &eval_grid)?,
      ConvolutionAlgorithm::Ndarray => ndarray_convolve(&f, &g, &eval_grid)?,
    };

    // Compute analytical expected result
    let expected_result = exponential_convolution(test_case.a, test_case.b, test_case.eval_domain, test_case.dx)?;

    // Compute metrics
    let metrics = DomainAgreementMetrics::new(actual_result.x(), actual_result.y(), expected_result.y())?;

    let execution_time = start_time.elapsed().as_secs_f64() * 1000.0;

    // Find peak error location
    let abs_errors: Array1<f64> = actual_result
      .y()
      .iter()
      .zip(expected_result.y().iter())
      .map(|(&a, &e)| (a - e).abs())
      .collect();

    let max_error_idx = abs_errors
      .iter()
      .enumerate()
      .max_by_key(|&(_, &a)| OrderedFloat(a))
      .unwrap()
      .0;

    let peak_error_location = actual_result.x()[max_error_idx];
    let peak_error_value = abs_errors[max_error_idx];

    Ok(TestResult {
      test_case_name: test_case.name.clone(),
      algorithm,
      test_case: test_case.clone(),
      metrics,
      execution_time_ms: execution_time,
      grid_points: eval_grid.len(),
      peak_error_location,
      peak_error_value,
    })
  }

  fn test_cases(&self) -> &[ExponentialTestCase] {
    &self.test_cases
  }

  fn function_type(&self) -> &'static str {
    "exponential"
  }
}

/// Exponential-specific flat result for TSV output
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExponentialFlatResult {
  pub test_case_name: String,
  pub algorithm: String,
  pub a: f64,
  pub b: f64,
  pub dx: f64,
  pub grid_points: usize,
  pub execution_time_ms: f64,
  pub r_squared: f64,
  pub max_abs_error: f64,
  pub mean_abs_error: f64,
  pub max_rel_error_percent: f64,
  pub mean_rel_error_percent: f64,
  pub rmse: f64,
  pub correlation: f64,
  pub mass_conservation_error: f64,
  pub peak_error_location: f64,
  pub peak_error_value: f64,
  pub tolerance_strict_abs: f64,
  pub tolerance_moderate_abs: f64,
  pub tolerance_loose_abs: f64,
  pub tolerance_strict_rel: f64,
  pub tolerance_moderate_rel: f64,
  pub tolerance_loose_rel: f64,
  pub stress_type: String,
  pub overall_assessment: String,
}

impl ToFlatResult for TestResult<ExponentialTestCase> {
  type FlatResult = ExponentialFlatResult;

  fn to_flat_result(&self) -> Self::FlatResult {
    ExponentialFlatResult {
      test_case_name: self.test_case_name.clone(),
      algorithm: self.algorithm.to_string(),
      a: self.test_case.a,
      b: self.test_case.b,
      dx: self.test_case.dx,
      grid_points: self.grid_points,
      execution_time_ms: self.execution_time_ms,
      r_squared: self.metrics.quality_metrics.r_squared,
      max_abs_error: self.metrics.abs_error_stats.max,
      mean_abs_error: self.metrics.abs_error_stats.mean,
      max_rel_error_percent: self.metrics.rel_error_stats.max * 100.0,
      mean_rel_error_percent: self.metrics.rel_error_stats.mean * 100.0,
      rmse: self.metrics.quality_metrics.rmse,
      correlation: self.metrics.quality_metrics.correlation,
      mass_conservation_error: self.metrics.quality_metrics.mass_error,
      peak_error_location: self.peak_error_location,
      peak_error_value: self.peak_error_value,
      tolerance_strict_abs: self.metrics.abs_tolerance_fraction(0),
      tolerance_moderate_abs: self.metrics.abs_tolerance_fraction(1),
      tolerance_loose_abs: self.metrics.abs_tolerance_fraction(2),
      tolerance_strict_rel: self.metrics.rel_tolerance_fraction(0),
      tolerance_moderate_rel: self.metrics.rel_tolerance_fraction(1),
      tolerance_loose_rel: self.metrics.rel_tolerance_fraction(2),
      stress_type: self.test_case.stress_type.clone(),
      overall_assessment: self.metrics.overall_assessment().to_string(),
    }
  }
}
