use crate::distribution::reference::convolution_test::algorithm_summary::AlgorithmSummary;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestSummary {
  pub function_type: String,
  pub total_tests: usize,
  pub total_successes: usize,
  pub total_failures: usize,
  pub total_algorithms: usize,
  pub execution_time_total_ms: f64,
  pub algorithm_summaries: Vec<AlgorithmSummary>,
  pub overall_assessment: String,
}
