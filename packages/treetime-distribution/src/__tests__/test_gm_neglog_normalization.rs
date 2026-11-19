#[cfg(test)]
mod tests {
  use crate::Distribution;
  use crate::distribution_multiplication;
  use eyre::Report;
  use helpers::{GoldenInput, GoldenOutput};
  use ndarray::Array1;
  use std::collections::BTreeMap;
  use std::path::Path;
  use treetime_grid::MaxGridPoints;
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

    let child = Distribution::function(
      Array1::from_vec(input.time_points.clone()),
      Array1::from_vec(input.child_neglog.clone()),
    )?;
    let coalescent = Distribution::function(
      Array1::from_vec(input.time_points.clone()),
      Array1::from_vec(input.coalescent_neglog.clone()),
    )?;
    let actual = distribution_multiplication(&child, &coalescent, MaxGridPoints::default())?.normalize()?;

    pretty_assert_ulps_eq!(Array1::from_vec(expected.time_points.clone()), actual.t(), max_ulps = 4);
    pretty_assert_abs_diff_eq!(
      Array1::from_vec(expected.probabilities_relative.clone()).mapv(|p| -p.ln()),
      actual.y(),
      epsilon = 1e-11,
      "costs relative to the peak are -ln of the v0 probabilities relative to the peak; v0 rounds its underflow rescaling to 7.3e-12 nats",
    );
    pretty_assert_ulps_eq!(expected.peak_position, actual.likely_time()?.unwrap(), max_ulps = 4);

    Ok(())
  }

  mod helpers {
    use deser::Deserialize;

    #[derive(Deserialize)]
    pub(super) struct GoldenInput {
      pub(super) time_points: Vec<f64>,
      pub(super) child_neglog: Vec<f64>,
      pub(super) coalescent_neglog: Vec<f64>,
    }

    #[derive(Deserialize)]
    pub(super) struct GoldenOutput {
      pub(super) time_points: Vec<f64>,
      pub(super) probabilities_relative: Vec<f64>,
      pub(super) peak_position: f64,
    }
  }
}
