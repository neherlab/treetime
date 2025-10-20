use crate::distribution::reference::convolution_test::framework::results::TestResult;
use crate::distribution::reference::convolution_test::output::ToFlatResult;
use serde::{Deserialize, Serialize};

use super::test_cases::GaussianTestCase;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianFlatResult {
  pub test_case_name: String,
  pub algorithm: String,
  pub sigma_f: f64,
  pub sigma_g: f64,
  pub mu: f64,
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
  pub peak_error_x: f64,
  pub peak_error_abs: f64,
  pub tolerance_strict_abs: f64,
  pub tolerance_moderate_abs: f64,
  pub tolerance_loose_abs: f64,
  pub tolerance_strict_rel: f64,
  pub tolerance_moderate_rel: f64,
  pub tolerance_loose_rel: f64,
  pub stress_type: String,
  pub overall_assessment: String,
}

impl ToFlatResult for TestResult<GaussianTestCase> {
  type FlatResult = GaussianFlatResult;

  fn to_flat_result(&self) -> Self::FlatResult {
    GaussianFlatResult {
      test_case_name: self.test_case.name.clone(),
      algorithm: self.algorithm.to_string(),
      sigma_f: self.test_case.sigma_f,
      sigma_g: self.test_case.sigma_g,
      mu: self.test_case.mu,
      dx: self.test_case.dx,
      grid_points: self.evaluation_grid.len(),
      execution_time_ms: self.execution_time_ms,
      r_squared: self.metrics.aggregate.domain_agreement.quality_metrics.r_squared,
      max_abs_error: self.metrics.aggregate.domain_agreement.abs_error_stats.max,
      mean_abs_error: self.metrics.aggregate.domain_agreement.abs_error_stats.mean,
      max_rel_error_percent: self.metrics.aggregate.domain_agreement.rel_error_stats.max * 100.0,
      mean_rel_error_percent: self.metrics.aggregate.domain_agreement.rel_error_stats.mean * 100.0,
      rmse: self.metrics.aggregate.domain_agreement.quality_metrics.rmse,
      correlation: self.metrics.aggregate.domain_agreement.quality_metrics.correlation,
      mass_conservation_error: self.metrics.aggregate.domain_agreement.quality_metrics.mass_error,
      peak_error_x: self.metrics.aggregate.domain_agreement.max_error_location.x_value,
      peak_error_abs: self.metrics.aggregate.domain_agreement.abs_error_stats.max,
      tolerance_strict_abs: self.metrics.aggregate.domain_agreement.abs_tolerance_fraction(0),
      tolerance_moderate_abs: self.metrics.aggregate.domain_agreement.abs_tolerance_fraction(1),
      tolerance_loose_abs: self.metrics.aggregate.domain_agreement.abs_tolerance_fraction(2),
      tolerance_strict_rel: self.metrics.aggregate.domain_agreement.rel_tolerance_fraction(0),
      tolerance_moderate_rel: self.metrics.aggregate.domain_agreement.rel_tolerance_fraction(1),
      tolerance_loose_rel: self.metrics.aggregate.domain_agreement.rel_tolerance_fraction(2),
      stress_type: self.test_case.stress_type.clone(),
      overall_assessment: self.metrics.aggregate.domain_agreement.overall_assessment().to_string(),
    }
  }
}
