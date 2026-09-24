#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionNegLog;
  use crate::distribution_multiplication;
  use eyre::Report;
  use helpers::neglog_to_plain_normalized;
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
    let actual = neglog_to_plain_normalized(&distribution_multiplication(&child, &coalescent)?);

    pretty_assert_ulps_eq!(Array1::from_vec(expected.time_points.clone()), actual.t(), max_ulps = 4);
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

  mod helpers {
    use crate::__tests__::aliases::{DistributionNegLog, DistributionPlain};
    use crate::Distribution;
    use crate::distribution_core::function::DistributionFunction;
    use ndarray_stats::QuantileExt;

    pub(super) fn neglog_to_plain_normalized(distribution: &DistributionNegLog) -> DistributionPlain {
      let Distribution::Function(function) = distribution else {
        panic!("expected a function distribution, got {distribution:?}");
      };
      let minimum = *function.y().min().unwrap();
      if !minimum.is_finite() {
        return Distribution::Empty;
      }
      let values = function.y().mapv(|value| (minimum - value).exp());
      Distribution::Function(
        DistributionFunction::from_start_dx_values(function.x_min(), function.dx(), values).unwrap(),
      )
    }

    mod tests {
      use super::neglog_to_plain_normalized;
      use crate::__tests__::aliases::{DistributionNegLog, DistributionPlain};
      use ndarray::array;
      use treetime_utils::pretty_assert_ulps_eq;

      #[test]
      fn test_gm_neglog_normalization_helper_preserves_likelihood_ratios() {
        let distribution = DistributionNegLog::function(array![0.0, 1.0, 2.0], array![1004.0, 1000.0, 1003.0]).unwrap();
        let actual = neglog_to_plain_normalized(&distribution);
        let expected = array![(-4.0_f64).exp(), 1.0, (-3.0_f64).exp()];
        pretty_assert_ulps_eq!(expected, actual.y(), max_ulps = 4);
      }

      #[test]
      fn test_gm_neglog_normalization_helper_rejects_nonfinite_minimum() {
        let distribution = DistributionNegLog::function(
          array![0.0, 1.0, 2.0],
          array![f64::INFINITY, f64::INFINITY, f64::INFINITY],
        )
        .unwrap();
        assert_eq!(DistributionPlain::Empty, neglog_to_plain_normalized(&distribution));
      }
    }
  }
}
