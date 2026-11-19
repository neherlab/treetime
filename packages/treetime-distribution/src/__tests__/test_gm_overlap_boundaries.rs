#[cfg(test)]
mod tests {
  use crate::{distribution_division, distribution_multiplication};
  use eyre::Report;
  use helpers::{DistributionKind, GoldenInput, GoldenOutput, Operation};
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::path::Path;
  use treetime_grid::MaxGridPoints;
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

    let left = input.left.to_distribution()?;
    let right = input.right.to_distribution()?;
    let actual = match input.operation {
      Operation::Divide => distribution_division(&left, &right, MaxGridPoints::default())?,
      Operation::Multiply => distribution_multiplication(&left, &right, MaxGridPoints::default())?,
    };

    pretty_assert_eq!(expected.kind, DistributionKind::from_distribution(&actual));
    let [expected_start, expected_end] = expected.bounds;
    let (actual_start, actual_end) = actual.time_bounds().unwrap();
    pretty_assert_ulps_eq!(expected_start, actual_start, max_ulps = 4);
    pretty_assert_ulps_eq!(expected_end, actual_end, max_ulps = 4);

    if let Some([expected_start_value, expected_end_value]) = expected.endpoint_values {
      let actual = [actual.eval(expected_start)?, actual.eval(expected_end)?];
      pretty_assert_abs_diff_eq!(expected_start_value, actual[0], epsilon = 1e-12);
      pretty_assert_abs_diff_eq!(expected_end_value, actual[1], epsilon = 1e-12);
    }

    Ok(())
  }

  mod helpers {
    use crate::Distribution;
    use deser::Deserialize;
    use eyre::Report;
    use ndarray::Array1;

    #[derive(Clone, Copy, Debug, Deserialize)]
    #[deser(rename_all = "lowercase")]
    pub(super) enum Operation {
      Divide,
      Multiply,
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
    #[deser(rename_all = "lowercase")]
    pub(super) enum DistributionKind {
      Function,
      Point,
      Range,
    }

    impl DistributionKind {
      pub(super) fn from_distribution(distribution: &Distribution) -> Self {
        match distribution {
          Distribution::Function(_) => Self::Function,
          Distribution::Point(_) => Self::Point,
          Distribution::Range(_) => Self::Range,
          Distribution::Empty => {
            unreachable!("golden overlap cases produce concrete non-empty distributions")
          },
        }
      }
    }

    #[derive(Deserialize)]
    pub(super) struct GoldenInput {
      pub(super) operation: Operation,
      pub(super) left: GoldenOperand,
      pub(super) right: GoldenOperand,
    }

    #[derive(Deserialize)]
    pub(super) struct GoldenOperand {
      kind: DistributionKind,
      x: Vec<f64>,
      y: Vec<f64>,
    }

    impl GoldenOperand {
      pub(super) fn to_distribution(&self) -> Result<Distribution, Report> {
        match self.kind {
          DistributionKind::Function => {
            Distribution::function(Array1::from_vec(self.x.clone()), Array1::from_vec(self.y.clone()))
          },
          DistributionKind::Range => Ok(Distribution::range((self.x[0], self.x[self.x.len() - 1]), self.y[0])),
          DistributionKind::Point => Ok(Distribution::point(self.x[0], self.y[0])),
        }
      }
    }

    #[derive(Deserialize)]
    pub(super) struct GoldenOutput {
      pub(super) kind: DistributionKind,
      pub(super) bounds: [f64; 2],
      pub(super) endpoint_values: Option<[f64; 2]>,
    }
  }
}
