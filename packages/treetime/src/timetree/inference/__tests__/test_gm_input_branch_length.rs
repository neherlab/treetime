#[cfg(test)]
mod tests {
  use crate::timetree::inference::input_branch_length::{
    InputBranchLikelihood, InputBranchModel, saturated_difference,
  };
  use deser::Deserialize;
  use eyre::Report;
  use itertools::{Itertools, izip};
  use ndarray::Array1;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::path::Path;
  use std::sync::LazyLock;
  use treetime_distribution::Distribution;
  use treetime_grid::{BoundaryBehavior, MaxGridPoints};
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_relative_eq, pretty_assert_ulps_eq};

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

  #[test]
  fn test_gm_input_branch_length_poisson_first_cell_of_a_high_cost_matches_v0_to_its_rounding() -> Result<(), Report> {
    let name = "poisson_below_threshold";
    let expected = &OUTPUTS[name];
    let production = helpers::production(name)?;
    let cell_end = production.function().dx() * helpers::BRANCH_RATE;
    let at_zero = production.distribution.eval(0.0)?;
    let (lengths, reference): (Vec<f64>, Vec<f64>) = izip!(&expected.branch_lengths, &expected.neg_log_lh)
      .filter(|&(&length, _)| length < cell_end)
      .map(|(&length, &neg_log_lh)| (length, neg_log_lh - expected.neg_log_lh[0] + at_zero))
      .unzip();
    let actual = lengths
      .iter()
      .map(|length| production.distribution.eval(length / helpers::BRANCH_RATE))
      .collect::<Result<Array1<f64>, Report>>()?;

    assert!(
      lengths.len() > 2,
      "the v0 grid of {name} has points inside the first cell"
    );
    pretty_assert_ulps_eq!(Array1::from_vec(reference), actual, max_ulps = 4);
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
  #[case::zero_length(             "zero_length",             Some(0.0))]
  #[case::poisson_one_mutation(    "poisson_one_mutation",    Some(1.0))]
  #[case::poisson_below_threshold( "poisson_below_threshold", Some(499.0))]
  #[case::gaussian_at_threshold(   "gaussian_at_threshold",   None)]
  #[case::gaussian_very_long(      "gaussian_very_long",      None)]
  #[trace]
  fn test_gm_input_branch_length_distribution_starts_hard_at_zero(
    #[case] name: &str,
    #[case] mutations: Option<f64>,
  ) -> Result<(), Report> {
    let input = &INPUTS[name];
    let production = helpers::production(name)?;
    let function = production.function();
    let one_mutation = input.sequence_length.to_f64_exact().recip();

    pretty_assert_abs_diff_eq!(0.0, function.x_min(), epsilon = 1e-15);
    match (mutations, function.left_extrap()) {
      (None, actual) => assert_eq!(BoundaryBehavior::Hard, actual),
      (Some(mutations), BoundaryBehavior::HardEdgeCell(law)) => {
        pretty_assert_abs_diff_eq!(0.0, law.t_edge, epsilon = 1e-15);
        pretty_assert_relative_eq!(1e-3 * one_mutation / helpers::BRANCH_RATE, law.offset, max_relative = 1e-15);
        pretty_assert_relative_eq!(mutations, law.b, epsilon = 1e-12, max_relative = 1e-12);
      },
      (Some(_), actual) => panic!("the Poisson case {name} must start with an edge cell, got {actual:?}"),
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::gaussian_at_threshold("gaussian_at_threshold")]
  #[case::gaussian_skewed_pi(   "gaussian_skewed_pi")]
  #[case::gaussian_very_long(   "gaussian_very_long")]
  #[trace]
  fn test_gm_input_branch_length_gaussian_differs_from_v0_within_the_linear_interpolation_bound(
    #[case] name: &str,
  ) -> Result<(), Report> {
    let production = helpers::production(name)?;
    let InputBranchLikelihood::Gaussian { variance, .. } = production.likelihood else {
      panic!("the case {name} has a Gaussian likelihood");
    };
    let step = production.function().dx() * helpers::BRANCH_RATE;
    let bound = step * step / (8.0 * variance);

    for (length, difference) in helpers::differences_from_v0(name, &production)? {
      assert!(
        difference <= bound,
        "the difference {difference:e} from v0 at the branch length {length:e} exceeds the bound {bound:e}"
      );
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::below_one_tenth_mutation("below_one_tenth_mutation")]
  #[case::poisson_one_mutation(    "poisson_one_mutation")]
  #[case::poisson_ten_mutations(   "poisson_ten_mutations")]
  #[case::poisson_long_genome(     "poisson_long_genome")]
  #[case::poisson_below_threshold( "poisson_below_threshold")]
  #[trace]
  fn test_gm_input_branch_length_poisson_differs_from_v0_within_the_bound_of_each_cell(
    #[case] name: &str,
  ) -> Result<(), Report> {
    let input = &INPUTS[name];
    let production = helpers::production(name)?;
    let step = production.function().dx() * helpers::BRANCH_RATE;
    let mutations = input.mutation_length * input.sequence_length.to_f64_exact();
    let grid_end = production.function().x_max() * helpers::BRANCH_RATE;

    for (length, difference) in helpers::differences_from_v0(name, &production)? {
      let cell = (length / step).floor();
      if length > grid_end || cell < 1.0 {
        continue;
      }
      let bound = mutations / (8.0 * cell * cell);
      assert!(
        difference <= bound,
        "the difference {difference:e} from v0 at the branch length {length:e} in cell {cell} exceeds the bound {bound:e}"
      );
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::below_one_tenth_mutation("below_one_tenth_mutation")]
  #[case::poisson_one_mutation(    "poisson_one_mutation")]
  #[case::poisson_ten_mutations(   "poisson_ten_mutations")]
  #[case::poisson_long_genome(     "poisson_long_genome")]
  #[trace]
  fn test_gm_input_branch_length_poisson_first_cell_matches_v0(#[case] name: &str) -> Result<(), Report> {
    let expected = &OUTPUTS[name];
    let production = helpers::production(name)?;
    let cell_end = production.function().dx() * helpers::BRANCH_RATE;
    let at_zero = production.distribution.eval(0.0)?;
    let (lengths, reference): (Vec<f64>, Vec<f64>) = izip!(&expected.branch_lengths, &expected.neg_log_lh)
      .filter(|&(&length, _)| length < cell_end)
      .map(|(&length, &neg_log_lh)| (length, neg_log_lh - expected.neg_log_lh[0]))
      .unzip();
    let actual = lengths
      .iter()
      .map(|length| Ok(production.distribution.eval(length / helpers::BRANCH_RATE)? - at_zero))
      .collect::<Result<Array1<f64>, Report>>()?;

    assert!(lengths.len() > 2, "the v0 grid of {name} has points inside the first cell");
    pretty_assert_abs_diff_eq!(Array1::from_vec(reference), actual, epsilon = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::gaussian_at_threshold("gaussian_at_threshold")]
  #[case::gaussian_skewed_pi(   "gaussian_skewed_pi")]
  #[case::gaussian_very_long(   "gaussian_very_long")]
  #[trace]
  fn test_gm_input_branch_length_gaussian_interpolation_error_converges_at_second_order(
    #[case] name: &str,
  ) -> Result<(), Report> {
    let input = &INPUTS[name];
    let likelihood = InputBranchModel::new(input.sequence_length, OUTPUTS[name].saturated_difference)
      .likelihood(input.mutation_length);
    let InputBranchLikelihood::Gaussian {
      mutation_length,
      variance,
      ..
    } = likelihood
    else {
      panic!("the case {name} has a Gaussian likelihood");
    };
    let quadratic = (mutation_length - 5.0 * variance.sqrt(), mutation_length + 5.0 * variance.sqrt());

    let normalized_errors = [10.0, 20.0, 40.0]
      .into_iter()
      .map(|points_per_width| {
        let distribution =
          likelihood.distribution_with_points_per_width(helpers::BRANCH_RATE, points_per_width, MaxGridPoints::default())?;
        let Distribution::Function(function) = distribution.as_ref() else {
          panic!("an input branch-length distribution is a Function");
        };
        let step = function.dx() * helpers::BRANCH_RATE;
        let at_zero = distribution.eval(0.0)? - likelihood.neg_log_lh(0.0);
        let largest = helpers::sub_cell_lengths(function, quadratic)
          .into_iter()
          .map(|length| {
            let actual = distribution.eval(length / helpers::BRANCH_RATE)? - at_zero;
            Ok::<f64, Report>((actual - likelihood.neg_log_lh(length)).abs())
          })
          .process_results(|errors| errors.fold(0.0, f64::max))?;
        Ok(largest * 8.0 * variance / (step * step))
      })
      .collect::<Result<Array1<f64>, Report>>()?;

    pretty_assert_relative_eq!(Array1::<f64>::ones(3), normalized_errors, max_relative = 1e-6);
    Ok(())
  }

  #[test]
  #[cfg_attr(
    dylint_lib = "custom",
    expect(
      versioned_name,
      reason = "`v0` names TreeTime 0.x, the reference implementation of this golden master, not a revision"
    )
  )]
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
    use treetime_distribution::{DistributionFunction, NegLog};
    use treetime_utils::cast::f64_to_usize_saturating;

    pub(super) const BRANCH_RATE: f64 = 0.002;

    #[cfg_attr(
      dylint_lib = "custom",
      expect(
        versioned_name,
        reason = "`v0` names TreeTime 0.x, the reference implementation of this golden master, not a revision"
      )
    )]
    pub(super) fn differences_from_v0(name: &str, production: &Production) -> Result<Vec<(f64, f64)>, Report> {
      let expected = &OUTPUTS[name];
      let at_zero = production.distribution.eval(0.0)?;
      izip!(&expected.branch_lengths, &expected.neg_log_lh)
        .map(|(&length, &neg_log_lh)| {
          let actual = production.distribution.eval(length / BRANCH_RATE)?;
          Ok((length, (actual - (neg_log_lh - expected.neg_log_lh[0] + at_zero)).abs()))
        })
        .collect()
    }

    pub(super) fn sub_cell_lengths(function: &DistributionFunction<f64, NegLog>, (start, end): (f64, f64)) -> Vec<f64> {
      let step = function.dx() * BRANCH_RATE;
      let first = f64_to_usize_saturating((start / step).ceil());
      let last = f64_to_usize_saturating((end / step).floor());
      (first..last)
        .flat_map(|cell| (0..64).map(move |k| (cell.to_f64_exact() + k.to_f64_exact() / 64.0) * step))
        .collect()
    }

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
