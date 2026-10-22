#[cfg(test)]
mod tests {
  use crate::timetree::inference::input_branch_length::{
    InputBranchLikelihood, InputBranchModel, saturated_difference,
  };
  use deser::Deserialize;
  use eyre::Report;
  use ndarray::Array1;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::path::Path;
  use std::sync::LazyLock;
  use treetime_grid::BoundaryBehavior;
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_relative_eq};

  const FIXTURES_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/src/timetree/inference/__tests__/__fixtures__"
  );

  static INPUTS: LazyLock<BTreeMap<String, helpers::GoldenInput>> = LazyLock::new(|| {
    json_read_file(Path::new(FIXTURES_DIR).join("gm_input_branch_length_inputs.json"))
      .expect("golden-master inputs must parse")
  });

  static OUTPUTS: LazyLock<BTreeMap<String, helpers::GoldenOutput>> = LazyLock::new(|| {
    json_read_file(Path::new(FIXTURES_DIR).join("gm_input_branch_length_outputs.json"))
      .expect("golden-master outputs must parse")
  });

  #[rustfmt::skip]
  #[rstest]
  #[case::zero_length(             "zero_length")]
  #[case::below_one_tenth_mutation("below_one_tenth_mutation")]
  #[case::poisson_one_mutation(    "poisson_one_mutation")]
  #[case::poisson_ten_mutations(   "poisson_ten_mutations")]
  #[case::poisson_long_genome(     "poisson_long_genome")]
  #[case::poisson_below_threshold( "poisson_below_threshold")]
  #[case::gaussian_at_threshold(   "gaussian_at_threshold")]
  #[case::gaussian_skewed_pi(      "gaussian_skewed_pi")]
  #[case::gaussian_very_long(      "gaussian_very_long")]
  #[trace]
  fn test_gm_input_branch_length_neg_log_lh_matches_v0(#[case] name: &str) -> Result<(), Report> {
    let input = &INPUTS[name];
    let expected = &OUTPUTS[name];
    let difference = saturated_difference(&Array1::from_vec(input.pi.clone()));
    pretty_assert_relative_eq!(expected.saturated_difference, difference, max_relative = 1e-15);
    let model = InputBranchModel::new(input.sequence_length, difference);
    let likelihood = model.likelihood(input.mutation_length);

    let neg_log_lh = Array1::from_vec(expected.branch_lengths.clone()).mapv(|length| likelihood.neg_log_lh(length));
    let peak = neg_log_lh.fold(f64::INFINITY, |min, &value| min.min(value));
    let actual = neg_log_lh.mapv(|value| value - peak);

    pretty_assert_relative_eq!(
      Array1::from_vec(expected.neg_log_lh.clone()),
      actual,
      epsilon = 1e-9,
      max_relative = 1e-12,
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero_length(             "zero_length")]
  #[case::below_one_tenth_mutation("below_one_tenth_mutation")]
  #[case::poisson_one_mutation(    "poisson_one_mutation")]
  #[case::poisson_ten_mutations(   "poisson_ten_mutations")]
  #[case::poisson_long_genome(     "poisson_long_genome")]
  #[case::poisson_below_threshold( "poisson_below_threshold")]
  #[case::gaussian_at_threshold(   "gaussian_at_threshold")]
  #[case::gaussian_skewed_pi(      "gaussian_skewed_pi")]
  #[case::gaussian_very_long(      "gaussian_very_long")]
  #[trace]
  fn test_gm_input_branch_length_distribution_equals_v0_formula_at_grid_nodes(#[case] name: &str) -> Result<(), Report> {
    let production = helpers::production(name)?;
    let function = production.function();
    let lengths = Array1::linspace(function.x_min(), function.x_max(), function.values().len())
      .mapv(|t| t * helpers::BRANCH_RATE);

    let formula = lengths.mapv(|length| production.likelihood.neg_log_lh(length));
    let peak = formula.fold(f64::INFINITY, |min, &value| min.min(value));
    let expected = formula.mapv(|value| value - peak);
    let actual = lengths
      .iter()
      .map(|length| production.distribution.eval(length / helpers::BRANCH_RATE))
      .collect::<Result<Array1<f64>, Report>>()?;

    pretty_assert_relative_eq!(expected, actual, epsilon = 1e-9, max_relative = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::gaussian_at_threshold("gaussian_at_threshold")]
  #[case::gaussian_skewed_pi(   "gaussian_skewed_pi")]
  #[case::gaussian_very_long(   "gaussian_very_long")]
  #[trace]
  fn test_gm_input_branch_length_gaussian_grid_ends_where_the_v0_cost_turns_linear(
    #[case] name: &str,
  ) -> Result<(), Report> {
    let production = helpers::production(name)?;
    let InputBranchLikelihood::Gaussian {
      mutation_length,
      variance,
      tail_scale,
    } = production.likelihood
    else {
      panic!("the case {name} has a Gaussian likelihood");
    };
    let deviation = production.function().x_max() * helpers::BRANCH_RATE - mutation_length;

    pretty_assert_relative_eq!(
      0.5 * deviation * deviation / variance,
      100.0 + deviation / tail_scale,
      max_relative = 1e-12
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::gaussian_at_threshold("gaussian_at_threshold")]
  #[case::gaussian_skewed_pi(   "gaussian_skewed_pi")]
  #[case::gaussian_very_long(   "gaussian_very_long")]
  #[trace]
  fn test_gm_input_branch_length_gaussian_tail_has_the_v0_slope(#[case] name: &str) -> Result<(), Report> {
    let expected = &OUTPUTS[name];
    let production = helpers::production(name)?;
    let BoundaryBehavior::Linear(law) = production.function().right_extrap() else {
      panic!("the right tail of {name} must be linear");
    };

    let n = expected.branch_lengths.len();
    let v0_slope = (expected.neg_log_lh[n - 1] - expected.neg_log_lh[n - 2])
      / (expected.branch_lengths[n - 1] - expected.branch_lengths[n - 2]);
    pretty_assert_relative_eq!(v0_slope * helpers::BRANCH_RATE, law.slope, max_relative = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero_length(             "zero_length")]
  #[case::poisson_one_mutation(    "poisson_one_mutation")]
  #[case::poisson_below_threshold( "poisson_below_threshold")]
  #[case::gaussian_at_threshold(   "gaussian_at_threshold")]
  #[case::gaussian_very_long(      "gaussian_very_long")]
  #[trace]
  fn test_gm_input_branch_length_distribution_starts_hard_at_zero(#[case] name: &str) -> Result<(), Report> {
    let production = helpers::production(name)?;
    let function = production.function();

    pretty_assert_abs_diff_eq!(0.0, function.x_min(), epsilon = 1e-15);
    assert_eq!(BoundaryBehavior::Hard, function.left_extrap());
    Ok(())
  }

  #[test]
  fn test_gm_input_branch_length_zero_length_distribution_matches_v0() -> Result<(), Report> {
    let expected = &OUTPUTS["zero_length"];
    let production = helpers::production("zero_length")?;
    let reference = Array1::from_vec(expected.neg_log_lh.clone());
    let actual = expected
      .branch_lengths
      .iter()
      .map(|length| production.distribution.eval(length / helpers::BRANCH_RATE))
      .collect::<Result<Array1<f64>, Report>>()?;

    pretty_assert_relative_eq!(reference, actual, epsilon = 1e-9, max_relative = 1e-12);
    Ok(())
  }

  mod helpers {
    use super::*;
    use std::sync::Arc;
    use treetime_distribution::{Distribution, DistributionFunction, NegLog};
    use treetime_grid::MaxGridPoints;

    pub(super) const BRANCH_RATE: f64 = 0.002;

    pub(super) fn production(name: &str) -> Result<Production, Report> {
      let input = &INPUTS[name];
      let likelihood = InputBranchModel::new(input.sequence_length, OUTPUTS[name].saturated_difference)
        .likelihood(input.mutation_length);
      let distribution = likelihood.distribution(BRANCH_RATE, MaxGridPoints::default())?;
      Ok(Production {
        likelihood,
        distribution,
      })
    }

    #[derive(Deserialize)]
    pub(super) struct GoldenInput {
      pub(super) mutation_length: f64,
      pub(super) sequence_length: usize,
      pub(super) pi: Vec<f64>,
    }

    #[derive(Deserialize)]
    pub(super) struct GoldenOutput {
      pub(super) saturated_difference: f64,
      pub(super) branch_lengths: Vec<f64>,
      pub(super) neg_log_lh: Vec<f64>,
    }

    pub(super) struct Production {
      pub(super) likelihood: InputBranchLikelihood,
      pub(super) distribution: Arc<Distribution<NegLog>>,
    }

    impl Production {
      pub(super) fn function(&self) -> &DistributionFunction<f64, NegLog> {
        let Distribution::Function(function) = self.distribution.as_ref() else {
          panic!(
            "an input branch-length distribution must be a Function, got {:?}",
            self.distribution
          );
        };
        function
      }
    }
  }
}
