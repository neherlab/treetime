use crate::testing::metrics::config::PointwiseConfig;
use ndarray::Array1;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use treetime_utils::array::serde::{array1_as_vec, array1_from_vec};

fn serialize_array_of_array1<S>(arrays: &[Array1<f64>; 3], serializer: S) -> Result<S::Ok, S::Error>
where
  S: Serializer,
{
  use serde::ser::SerializeSeq;
  let mut seq = serializer.serialize_seq(Some(3))?;
  for array in arrays {
    let vec: Vec<f64> = array.to_vec();
    seq.serialize_element(&vec)?;
  }
  seq.end()
}

fn deserialize_array_of_array1<'de, D>(deserializer: D) -> Result<[Array1<f64>; 3], D::Error>
where
  D: Deserializer<'de>,
{
  let vecs: Vec<Vec<f64>> = Vec::deserialize(deserializer)?;
  if vecs.len() != 3 {
    return Err(serde::de::Error::custom(format!(
      "Expected 3 arrays, got {}",
      vecs.len()
    )));
  }
  Ok([
    Array1::from_vec(vecs[0].clone()),
    Array1::from_vec(vecs[1].clone()),
    Array1::from_vec(vecs[2].clone()),
  ])
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToleranceMetrics {
  #[serde(
    serialize_with = "serialize_array_of_array1",
    deserialize_with = "deserialize_array_of_array1"
  )]
  pub pass_masks: [Array1<f64>; 3],
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub support_coverage_mask: Array1<f64>,
  pub summary: ToleranceSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToleranceSummary {
  pub pass_fractions: [f64; 3],
  pub support_mismatch_count: usize,
}

pub(super) fn compute_tolerance_metrics(
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  config: &PointwiseConfig,
) -> eyre::Result<ToleranceMetrics> {
  let n = actual.len();

  let abs_errors = (actual - expected).mapv(|x| x.abs());
  let mut rel_errors = Array1::zeros(n);
  for i in 0..n {
    let denom = expected[i].abs().max(config.epsilon);
    rel_errors[i] = abs_errors[i] / denom;
  }

  let mut pass_masks = [Array1::zeros(n), Array1::zeros(n), Array1::zeros(n)];
  let mut pass_counts = [0_usize; 3];

  for level in 0..3 {
    let abs_tol = config.abs_tolerances[level];
    let rel_tol = config.rel_tolerances[level];
    for i in 0..n {
      let abs_pass = abs_errors[i] <= abs_tol;
      let rel_pass = rel_errors[i] <= rel_tol;
      if abs_pass || rel_pass {
        pass_masks[level][i] = 1.0;
        pass_counts[level] += 1;
      }
    }
  }

  let mut support_coverage_mask = Array1::zeros(n);
  let mut support_mismatch_count = 0;
  for i in 0..n {
    let actual_significant = actual[i].abs() > config.log_threshold;
    let expected_significant = expected[i].abs() > config.log_threshold;
    if actual_significant != expected_significant {
      support_coverage_mask[i] = 1.0;
      support_mismatch_count += 1;
    }
  }

  let pass_fractions = [
    pass_counts[0] as f64 / n as f64,
    pass_counts[1] as f64 / n as f64,
    pass_counts[2] as f64 / n as f64,
  ];

  let summary = ToleranceSummary {
    pass_fractions,
    support_mismatch_count,
  };

  Ok(ToleranceMetrics {
    pass_masks,
    support_coverage_mask,
    summary,
  })
}
