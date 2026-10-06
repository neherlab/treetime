use crate::testing::framework::test_case::TestCase;
use crate::testing::metrics::metrics::ValidationMetrics;
use deser::{Deserialize, Serialize};
use ndarray::Array1;
use treetime_utils::adapters::ArrayVec;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum TestRunOutcome<T: TestCase> {
  Success(Box<TestResult<T>>),
  Failure(TestFailure<T>),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestResult<T: TestCase> {
  pub(crate) algorithm: String,
  pub(crate) test_case: T,
  pub(crate) execution_time_ms: f64,

  #[deser(as = ArrayVec)]
  pub(crate) f_x_values: Array1<f64>,
  #[deser(as = ArrayVec)]
  pub(crate) f_y_values: Array1<f64>,
  #[deser(as = ArrayVec)]
  pub(crate) g_x_values: Array1<f64>,
  #[deser(as = ArrayVec)]
  pub(crate) g_y_values: Array1<f64>,

  #[deser(as = ArrayVec)]
  pub(crate) evaluation_grid: Array1<f64>,
  #[deser(as = ArrayVec)]
  pub(crate) actual_values: Array1<f64>,
  #[deser(as = ArrayVec)]
  pub(crate) expected_values: Array1<f64>,

  pub(crate) metrics: ValidationMetrics,

  #[deser(default, skip_serializing_if = Option::is_none)]
  pub(crate) log_scale_actual: Option<f64>,
  #[deser(default, skip_serializing_if = Option::is_none)]
  pub(crate) log_scale_expected: Option<f64>,
  #[deser(default, skip_serializing_if = Option::is_none)]
  pub(crate) log_scale_error: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestFailure<T: TestCase> {
  pub(crate) algorithm: String,
  pub(crate) test_case: T,
  pub(crate) error: String,
  pub(crate) execution_time_ms: f64,
}
