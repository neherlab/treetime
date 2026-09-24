use crate::testing::framework::test_case::TestCase;
use crate::testing::metrics::metrics::ValidationMetrics;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use treetime_utils::array::serde::{array1_as_vec, array1_from_vec};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum TestRunOutcome<T: TestCase> {
  Success(Box<TestResult<T>>),
  Failure(TestFailure<T>),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestResult<T: TestCase> {
  pub(crate) algorithm: String,
  pub(crate) test_case: T,
  pub(crate) execution_time_ms: f64,

  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub(crate) f_x_values: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub(crate) f_y_values: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub(crate) g_x_values: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub(crate) g_y_values: Array1<f64>,

  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub(crate) evaluation_grid: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub(crate) actual_values: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub(crate) expected_values: Array1<f64>,

  pub(crate) metrics: ValidationMetrics,

  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub(crate) log_scale_actual: Option<f64>,
  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub(crate) log_scale_expected: Option<f64>,
  #[serde(default, skip_serializing_if = "Option::is_none")]
  pub(crate) log_scale_error: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestFailure<T: TestCase> {
  pub(crate) algorithm: String,
  pub(crate) test_case: T,
  pub(crate) error: String,
  pub(crate) execution_time_ms: f64,
}
