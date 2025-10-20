use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::convolution_test::metrics::metrics::ConvolutionMetrics;
use crate::io::serde::{array1_as_vec, array1_from_vec};
use ndarray::Array1;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestResult<T: TestCase> {
  pub algorithm: String,
  pub test_case: T,
  pub execution_time_ms: f64,

  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub f_x_values: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub f_y_values: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub g_x_values: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub g_y_values: Array1<f64>,

  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub evaluation_grid: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub actual_values: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub expected_values: Array1<f64>,

  pub metrics: ConvolutionMetrics,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestFailure<T: TestCase> {
  pub algorithm: String,
  pub test_case: T,
  pub error: String,
  pub execution_time_ms: f64,
}

#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TestRunOutcome<T: TestCase> {
  Success(TestResult<T>),
  Failure(TestFailure<T>),
}
