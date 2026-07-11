#[cfg(test)]
mod tests {
  use crate::{DistributionNegLog, distribution_multiplication};
  use eyre::Report;
  use ndarray::Array1;
  use serde::Deserialize;
  use std::collections::BTreeMap;
  use std::path::Path;
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/src/__tests__/__fixtures__");

  #[test]
  fn test_gm_neglog_normalization_matches_v0_coalescent_underflow() -> Result<(), Report> {
    let inputs: BTreeMap<String, GoldenInput> =
      json_read_file(Path::new(FIXTURES_DIR).join("gm_neglog_normalization_inputs.json"))?;
    let outputs: BTreeMap<String, GoldenOutput> =
      json_read_file(Path::new(FIXTURES_DIR).join("gm_neglog_normalization_outputs.json"))?;
    let input = &inputs["coalescent_underflow"];
    let expected = &outputs["coalescent_underflow"];

    let child = DistributionNegLog::function(
      Array1::from_vec(input.time_points.clone()),
      Array1::from_vec(input.child_neglog.clone()),
    )?;
    let coalescent = DistributionNegLog::function(
      Array1::from_vec(input.time_points.clone()),
      Array1::from_vec(input.coalescent_neglog.clone()),
    )?;
    let actual = distribution_multiplication(&child, &coalescent)?.to_plain_normalized();

    pretty_assert_ulps_eq!(Array1::from_vec(expected.time_points.clone()), actual.t(), max_ulps = 4);
    // v0 clips endpoint evaluation inward by TINY_NUMBER=1e-12 before
    // interpolation (`treetime/distribution.py:295-303`).
    pretty_assert_abs_diff_eq!(
      Array1::from_vec(expected.probabilities_relative.clone()),
      actual.y(),
      epsilon = 1e-12,
    );
    pretty_assert_ulps_eq!(expected.peak_position, actual.likely_time().unwrap(), max_ulps = 4);

    Ok(())
  }

  #[derive(Deserialize)]
  struct GoldenInput {
    time_points: Vec<f64>,
    child_neglog: Vec<f64>,
    coalescent_neglog: Vec<f64>,
  }

  #[derive(Deserialize)]
  struct GoldenOutput {
    time_points: Vec<f64>,
    probabilities_relative: Vec<f64>,
    peak_position: f64,
  }
}
