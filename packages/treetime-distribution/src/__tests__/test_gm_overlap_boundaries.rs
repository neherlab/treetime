#[cfg(test)]
mod tests {
  use crate::{distribution_division, distribution_multiplication};
  use eyre::Report;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::path::Path;
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_eq, pretty_assert_ulps_eq};

  const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/src/__tests__/__fixtures__");

  #[rustfmt::skip]
  #[rstest]
  #[case::divide_endpoint_contact("divide_endpoint_contact")]
  #[case::divide_function_narrower_divisor("divide_function_narrower_divisor")]
  #[case::divide_range_partial_overlap("divide_range_partial_overlap")]
  #[case::multiply_endpoint_contact("multiply_endpoint_contact")]
  #[case::multiply_range_partial_overlap("multiply_range_partial_overlap")]
  #[trace]
  fn test_gm_overlap_boundaries_matches_v0(#[case] case_name: &str) -> Result<(), Report> {
    let inputs: BTreeMap<String, GoldenInput> =
      json_read_file(Path::new(FIXTURES_DIR).join("gm_overlap_boundaries_inputs.json"))?;
    let outputs: BTreeMap<String, GoldenOutput> =
      json_read_file(Path::new(FIXTURES_DIR).join("gm_overlap_boundaries_outputs.json"))?;
    let input = &inputs[case_name];
    let expected = &outputs[case_name];

    let left = input.left().to_distribution()?;
    let right = input.right().to_distribution()?;
    let actual = match input.operation() {
      Operation::Divide => distribution_division(&left, &right)?,
      Operation::Multiply => distribution_multiplication(&left, &right)?,
    };

    // Oracle: v0 `Distribution.multiply()` and `Distribution.divide()` in
    // packages/legacy/treetime/treetime/distribution.py:82-185.
    pretty_assert_eq!(expected.kind(), DistributionKind::from_distribution(&actual));
    let [expected_start, expected_end] = expected.bounds();
    let (actual_start, actual_end) = actual.time_bounds().unwrap();
    pretty_assert_ulps_eq!(expected_start, actual_start, max_ulps = 4);
    pretty_assert_ulps_eq!(expected_end, actual_end, max_ulps = 4);

    if let Some([expected_start_value, expected_end_value]) = expected.endpoint_values() {
      let actual = [actual.eval(expected_start)?, actual.eval(expected_end)?];
      // V0 clips vectorized endpoint evaluation inward by TINY_NUMBER=1e-12.
      pretty_assert_abs_diff_eq!(expected_start_value, actual[0], epsilon = 1e-12);
      pretty_assert_abs_diff_eq!(expected_end_value, actual[1], epsilon = 1e-12);
    }

    Ok(())
  }

  mod helpers {
    use crate::DistributionNegLog;
    use eyre::Report;
    use ndarray::Array1;
    use serde::Deserialize;

    #[derive(Clone, Copy, Debug, Deserialize)]
    #[serde(rename_all = "lowercase")]
    pub(super) enum Operation {
      Divide,
      Multiply,
    }

    #[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
    #[serde(rename_all = "lowercase")]
    pub(super) enum DistributionKind {
      Function,
      Point,
      Range,
    }

    impl DistributionKind {
      pub(super) fn from_distribution(distribution: &DistributionNegLog) -> Self {
        match distribution {
          DistributionNegLog::Function(_) => Self::Function,
          DistributionNegLog::Point(_) => Self::Point,
          DistributionNegLog::Range(_) => Self::Range,
          DistributionNegLog::Empty | DistributionNegLog::Formula(_) => {
            unreachable!("golden overlap cases produce concrete non-empty distributions")
          },
        }
      }
    }

    #[derive(Deserialize)]
    pub(super) struct GoldenInput {
      operation: Operation,
      left: GoldenOperand,
      right: GoldenOperand,
    }

    impl GoldenInput {
      pub(super) fn operation(&self) -> Operation {
        self.operation
      }

      pub(super) fn left(&self) -> &GoldenOperand {
        &self.left
      }

      pub(super) fn right(&self) -> &GoldenOperand {
        &self.right
      }
    }

    #[derive(Deserialize)]
    pub(super) struct GoldenOperand {
      kind: DistributionKind,
      x: Vec<f64>,
      y: Vec<f64>,
    }

    impl GoldenOperand {
      pub(super) fn to_distribution(&self) -> Result<DistributionNegLog, Report> {
        match self.kind {
          DistributionKind::Function => {
            DistributionNegLog::function(Array1::from_vec(self.x.clone()), Array1::from_vec(self.y.clone()))
          },
          DistributionKind::Range => Ok(DistributionNegLog::range(
            (self.x[0], self.x[self.x.len() - 1]),
            self.y[0],
          )),
          DistributionKind::Point => Ok(DistributionNegLog::point(self.x[0], self.y[0])),
        }
      }
    }

    #[derive(Deserialize)]
    pub(super) struct GoldenOutput {
      kind: DistributionKind,
      bounds: [f64; 2],
      endpoint_values: Option<[f64; 2]>,
    }

    impl GoldenOutput {
      pub(super) fn kind(&self) -> DistributionKind {
        self.kind
      }

      pub(super) fn bounds(&self) -> [f64; 2] {
        self.bounds
      }

      pub(super) fn endpoint_values(&self) -> Option<[f64; 2]> {
        self.endpoint_values
      }
    }
  }

  use helpers::{DistributionKind, GoldenInput, GoldenOutput, Operation};
}
