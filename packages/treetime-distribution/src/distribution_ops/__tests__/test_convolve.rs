#[cfg(test)]
mod tests {
  use self::helpers::{
    DistributionVariant, PiecewiseLinear, distribution, exponential_box_density, function_on, neg_log,
  };
  use crate::__tests__::aliases::DistributionNegLog;
  use crate::distribution_core::function::DistributionFunction;
  use crate::distribution_ops::__tests__::test_mass_profile_convergence::tests::helpers::assert_second_order;
  use crate::distribution_ops::convolve::{coarsen_convolution, distribution_convolution};
  use approx::ulps_eq;
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::{Array1, array, s};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use statrs::distribution::{ContinuousCDF, Normal};
  use std::f64::consts::{LN_2, PI};
  use treetime_grid::grid::Grid;
  use treetime_grid::{BoundaryBehavior, HardApproachLaw, MaxGridPoints, SoftTailLaw};
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  #[test]
  fn test_convolve_neglog_gaussian_variances_add() -> Result<(), Report> {
    let gaussian = |sigma: f64| -> Result<DistributionNegLog, Report> {
      let t = Array1::linspace(-6.0, 6.0, 241);
      let y = t.mapv(|t| t * t / (2.0 * sigma * sigma));
      DistributionNegLog::function(t, y)
    };

    let result = distribution_convolution(&gaussian(1.0)?, &gaussian(1.0)?, MaxGridPoints::default())?;

    let y0 = result.eval(0.0)?;
    pretty_assert_abs_diff_eq!(0.0, result.likely_time()?.unwrap(), epsilon = 0.05);
    pretty_assert_abs_diff_eq!(0.25, result.eval(1.0)? - y0, epsilon = 1e-3);
    pretty_assert_abs_diff_eq!(1.00, result.eval(2.0)? - y0, epsilon = 1e-3);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::formula_empty(   DistributionVariant::Formula,  DistributionVariant::Empty,    "Cannot convolve formula with empty: operation not implemented")]
  #[case::formula_point(   DistributionVariant::Formula,  DistributionVariant::Point,    "Cannot convolve formula with point: operation not implemented")]
  #[case::formula_range(   DistributionVariant::Formula,  DistributionVariant::Range,    "Cannot convolve formula with range: operation not implemented")]
  #[case::formula_function(DistributionVariant::Formula,  DistributionVariant::Function, "Cannot convolve formula with function: operation not implemented")]
  #[case::formula_formula( DistributionVariant::Formula,  DistributionVariant::Formula,  "Cannot convolve formula with formula: operation not implemented")]
  #[case::empty_formula(   DistributionVariant::Empty,    DistributionVariant::Formula,  "Cannot convolve empty with formula: operation not implemented")]
  #[case::point_formula(   DistributionVariant::Point,    DistributionVariant::Formula,  "Cannot convolve point with formula: operation not implemented")]
  #[case::range_formula(   DistributionVariant::Range,    DistributionVariant::Formula,  "Cannot convolve range with formula: operation not implemented")]
  #[case::function_formula(DistributionVariant::Function, DistributionVariant::Formula,  "Cannot convolve function with formula: operation not implemented")]
  #[case::range_range(     DistributionVariant::Range,    DistributionVariant::Range,    "Cannot convolve range with range: operation not implemented")]
  #[trace]
  fn test_convolve_unsupported_pairs_return_errors(
    #[case] left: DistributionVariant,
    #[case] right: DistributionVariant,
    #[case] expected: &str,
  ) {
    assert_error!(distribution_convolution(&distribution(left), &distribution(right), MaxGridPoints::default()), expected);
  }

  #[test]
  fn test_convolution_empty() -> Result<(), Report> {
    let a = DistributionNegLog::empty();
    let b = DistributionNegLog::function(array![], array![])?;
    let actual = distribution_convolution(&a, &b, MaxGridPoints::default())?;
    assert_eq!(DistributionNegLog::empty(), actual);
    Ok(())
  }

  #[test]
  fn test_convolution_point_point_adds_neg_log_amplitudes() -> Result<(), Report> {
    let a = DistributionNegLog::point(2.0, 3.0);
    let b = DistributionNegLog::point(5.0, 4.0);
    let actual = distribution_convolution(&a, &b, MaxGridPoints::default())?;
    assert_eq!(DistributionNegLog::point(7.0, 7.0), actual);
    Ok(())
  }

  #[test]
  fn test_convolution_point_function_shifts_values_tail_laws_and_level() -> Result<(), Report> {
    let p = DistributionNegLog::point(3.0, 2.0);
    let f = DistributionFunction::from_start_dx_values(1.0, 1.0, array![3.0, 1.0, 2.0, 4.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.5, b: 2.0 }))?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 0.5 }))?
      .with_amplitude(1.5);

    let actual = distribution_convolution(&p, &DistributionNegLog::Function(f), MaxGridPoints::default())?;

    let expected = DistributionFunction::from_start_dx_values(4.0, 1.0, array![3.0, 1.0, 2.0, 4.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 3.5, b: 2.0 }))?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 0.5 }))?
      .with_amplitude(3.5);
    assert_eq!(DistributionNegLog::Function(expected), actual);
    Ok(())
  }

  #[test]
  fn test_convolution_point_range() -> Result<(), Report> {
    let p = DistributionNegLog::point(3.0, 2.0);
    let r = DistributionNegLog::range((1.0, 4.0), 1.5);
    let actual = distribution_convolution(&p, &r, MaxGridPoints::default())?;
    assert_eq!(DistributionNegLog::range((4.0, 7.0), 3.5), actual);
    Ok(())
  }

  #[test]
  fn test_convolution_range_point() -> Result<(), Report> {
    let r = DistributionNegLog::range((1.0, 4.0), 1.5);
    let p = DistributionNegLog::point(3.0, 2.0);
    let actual = distribution_convolution(&r, &p, MaxGridPoints::default())?;
    assert_eq!(DistributionNegLog::range((4.0, 7.0), 3.5), actual);
    Ok(())
  }

  #[test]
  fn test_convolution_range_function_covers_the_full_support() -> Result<(), Report> {
    let r = DistributionNegLog::range((2.0, 6.0), LN_2);
    let f = DistributionFunction::from_start_dx_values(0.0, 2.0, Array1::zeros(6))?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?;

    let actual = distribution_convolution(&r, &DistributionNegLog::Function(f), MaxGridPoints::default())?;

    let DistributionNegLog::Function(g) = &actual else {
      panic!("expected a Function result, got {actual:?}");
    };
    let overlap_length = |t: f64| ((t - 2.0).min(10.0) - (t - 6.0).max(0.0)).max(0.0);
    let rise_cell = 4.0 / 8.0;
    let t = Array1::linspace(2.0, 16.0, 1 + 14 * 2);
    pretty_assert_ulps_eq!(rise_cell, g.dx(), max_ulps = 0);
    pretty_assert_ulps_eq!(t, g.t(), max_ulps = 4);
    pretty_assert_abs_diff_eq!(
      t.mapv(overlap_length) / 2.0,
      g.values().mapv(|y| (-y).exp()),
      epsilon = 1e-15
    );
    assert_eq!(
      (BoundaryBehavior::Hard, BoundaryBehavior::Hard),
      (g.left_extrap(), g.right_extrap())
    );
    for t in [2.3, 15.7] {
      pretty_assert_abs_diff_eq!(overlap_length(t) / 2.0, (-actual.eval(t)?).exp(), epsilon = 1e-15);
    }
    Ok(())
  }

  #[test]
  fn test_convolution_range_exponential_matches_the_trapezoid_model() -> Result<(), Report> {
    let (lambda, range) = (1.5_f64, (1.0, 2.0));
    let f = function_on((0.0, 4.0), 17, |t| lambda * t - lambda.ln())?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: lambda }))?;
    let model = PiecewiseLinear::of(&f, lambda);

    let actual = distribution_convolution(
      &DistributionNegLog::range(range, 0.0),
      &DistributionNegLog::Function(f),
      MaxGridPoints::default(),
    )?;

    let DistributionNegLog::Function(g) = &actual else {
      panic!("expected a Function result, got {actual:?}");
    };
    let y = g.values();
    assert!(
      y[0].is_infinite(),
      "zero density at the start of the support, got {}",
      y[0]
    );
    let expected = g.t().slice(s![1..]).mapv(|t| -model.box_mass(range, t).ln());
    pretty_assert_abs_diff_eq!(expected, y.slice(s![1..]).to_owned(), epsilon = 1e-12);
    assert_eq!(
      (
        BoundaryBehavior::Hard,
        BoundaryBehavior::Linear(SoftTailLaw { slope: lambda })
      ),
      (g.left_extrap(), g.right_extrap())
    );
    pretty_assert_abs_diff_eq!(-model.box_mass(range, 9.0).ln(), actual.eval(9.0)?, epsilon = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::rising(  1.5)]
  #[case::falling( 2.5)]
  #[case::tail(    3.0)]
  #[trace]
  fn test_convolution_range_exponential_converges_at_second_order(#[case] t: f64) -> Result<(), Report> {
    let (lambda, range) = (1.5_f64, (1.0, 2.0));
    let errors: Vec<f64> = [16, 32, 64, 128]
      .iter()
      .map(|&cells| {
        let f = function_on((0.0, 4.0), cells + 1, |t| lambda * t - lambda.ln())?
          .with_left_extrap(BoundaryBehavior::Hard)?
          .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: lambda }))?;
        let actual = distribution_convolution(
          &DistributionNegLog::range(range, 0.0),
          &DistributionNegLog::Function(f),
          MaxGridPoints::default(),
        )?;
        helpers::assert_grid_point(&actual, t);
        Ok((-actual.eval(t)?).exp() - exponential_box_density(lambda, range, t))
      })
      .try_collect::<_, _, Report>()?;
    assert_second_order(&errors);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::rising(  0.5)]
  #[case::plateau( 1.0)]
  #[case::falling( 2.5)]
  #[case::tail(    4.0)]
  #[trace]
  fn test_convolution_range_gaussian_converges_at_second_order(#[case] t: f64) -> Result<(), Report> {
    let (a, b) = (0.0, 2.0);
    let normal = Normal::standard();
    let errors: Vec<f64> = [128, 256, 512, 1024]
      .iter()
      .map(|&cells| {
        let f = function_on((-8.0, 8.0), cells + 1, |x| x * x / 2.0 + (2.0 * PI).sqrt().ln())?
          .with_left_extrap(BoundaryBehavior::Hard)?
          .with_right_extrap(BoundaryBehavior::Hard)?;
        let actual = distribution_convolution(
          &DistributionNegLog::range((a, b), 0.0),
          &DistributionNegLog::Function(f),
          MaxGridPoints::default(),
        )?;
        helpers::assert_grid_point(&actual, t);
        Ok((-actual.eval(t)?).exp() - (normal.cdf(t - a) - normal.cdf(t - b)))
      })
      .try_collect::<_, _, Report>()?;
    assert_second_order(&errors);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::rising(  0.5)]
  #[case::peak(    1.0)]
  #[case::falling( 2.5)]
  #[trace]
  fn test_convolution_function_function_hard_edges_converge_at_second_order(#[case] t: f64) -> Result<(), Report> {
    let (a, b) = (1.5_f64, 4.0_f64);
    let exponential = |rate: f64, cells: usize| -> Result<DistributionNegLog, Report> {
      let f = function_on((0.0, 4.0), cells + 1, |x| rate * x - rate.ln())?
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: rate }))?;
      Ok(DistributionNegLog::Function(f))
    };
    let errors: Vec<f64> = [32, 64, 128, 256]
      .iter()
      .map(|&cells| {
        let actual = distribution_convolution(&exponential(a, cells)?, &exponential(b, cells)?, MaxGridPoints::default())?;
        helpers::assert_grid_point(&actual, t);
        let expected = a * b / (b - a) * ((-a * t).exp() - (-b * t).exp());
        Ok((-actual.eval(t)?).exp() - expected)
      })
      .try_collect::<_, _, Report>()?;
    assert_second_order(&errors);
    Ok(())
  }

  #[test]
  fn test_convolution_function_function_hard_edges_mass_converges_at_second_order() -> Result<(), Report> {
    let (a, b, end) = (1.5_f64, 4.0_f64, 4.0_f64);
    let exponential = |rate: f64, cells: usize| -> Result<DistributionNegLog, Report> {
      let f = function_on((0.0, end), cells + 1, |x| rate * x - rate.ln())?
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Hard)?;
      Ok(DistributionNegLog::Function(f))
    };
    let operand_mass = |rate: f64| 1.0 - (-rate * end).exp();
    let errors: Vec<f64> = [32, 64, 128, 256]
      .iter()
      .map(|&cells| {
        let actual = distribution_convolution(
          &exponential(a, cells)?,
          &exponential(b, cells)?,
          MaxGridPoints::default(),
        )?;
        let DistributionNegLog::Function(conv) = &actual else {
          panic!("expected a Function result, got {actual:?}");
        };
        let density = conv.values().mapv(|v| (-v).exp());
        let trapezoid_mass = conv.dx() * (density.sum() - 0.5 * (density[0] + density[density.len() - 1]));
        Ok(trapezoid_mass - operand_mass(a) * operand_mass(b))
      })
      .try_collect::<_, _, Report>()?;
    assert_second_order(&errors);
    Ok(())
  }

  #[test]
  fn test_convolution_range_function_carries_both_amplitudes() -> Result<(), Report> {
    let f = function_on((0.0, 3.0), 7, |t| t)?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 }))?;
    let model = PiecewiseLinear::of(&f, 1.0);
    let height = 3.0_f64;
    let f = f.with_amplitude(-2.0_f64.ln());

    let actual = distribution_convolution(
      &DistributionNegLog::range((0.0, 1.0), -height.ln()),
      &DistributionNegLog::Function(f),
      MaxGridPoints::default(),
    )?;

    let DistributionNegLog::Function(g) = &actual else {
      panic!("expected a Function result, got {actual:?}");
    };
    let expected = g
      .t()
      .slice(s![1..])
      .mapv(|t| -(2.0 * height * model.box_mass((0.0, 1.0), t)).ln());
    pretty_assert_abs_diff_eq!(expected, g.values().slice(s![1..]).to_owned(), epsilon = 1e-12);
    let beyond_grid = |t: f64| 6.0 * ((1.0 - t).exp() - (-t).exp());
    for t in [4.0, 5.0] {
      pretty_assert_abs_diff_eq!(-beyond_grid(t).ln(), actual.eval(t)?, epsilon = 1e-12);
    }
    Ok(())
  }

  #[test]
  fn test_convolution_range_function_keeps_far_tail_mass_against_large_amplitude() -> Result<(), Report> {
    let f = DistributionFunction::from_range_values((0.0, 1000.0), array![0.0, 1000.0])?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 }))?;
    let r = DistributionNegLog::range((0.0, 1.0), -1000.0);

    let actual = distribution_convolution(&r, &DistributionNegLog::Function(f), MaxGridPoints::default())?;

    let tail_cell_neg_log_mass = -(-(-1.0_f64).exp_m1()).ln();
    let cancelled_at_ulp_of_1000 = actual.eval(1001.0)?;
    pretty_assert_abs_diff_eq!(tail_cell_neg_log_mass, cancelled_at_ulp_of_1000, epsilon = 1e-13);
    Ok(())
  }

  #[test]
  fn test_convolution_range_function_box_narrower_than_spacing() -> Result<(), Report> {
    let f = DistributionFunction::from_range_values((0.0, 3.0), array![0.0, 1.0, 2.0, 3.0])?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 }))?;
    let model = PiecewiseLinear::of(&f, 1.0);

    let actual = distribution_convolution(
      &DistributionNegLog::range((0.0, 0.2), 0.0),
      &DistributionNegLog::Function(f),
      MaxGridPoints::default(),
    )?;

    let DistributionNegLog::Function(g) = &actual else {
      panic!("expected a Function result, got {actual:?}");
    };
    pretty_assert_ulps_eq!(Array1::linspace(0.0, 3.2, 1 + 128), g.t(), max_ulps = 4);
    let expected = g.t().slice(s![1..]).mapv(|t| -model.box_mass((0.0, 0.2), t).ln());
    pretty_assert_abs_diff_eq!(expected, g.y().slice(s![1..]).to_owned(), epsilon = 1e-15);
    assert!(g.y()[0].is_infinite(), "zero density at the start, got {}", g.y()[0]);
    Ok(())
  }

  #[test]
  fn test_convolution_range_function_end_cells_rise_linearly_from_the_support_ends() -> Result<(), Report> {
    let f = DistributionFunction::from_range_values((0.0, 1.0), array![0.0, 0.0])?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?;
    let box_mass = |t: f64| (t.min(1.0) - (t - 0.2).max(0.0)).max(0.0);

    let actual = distribution_convolution(
      &DistributionNegLog::range((0.0, 0.2), 0.0),
      &DistributionNegLog::Function(f),
      MaxGridPoints::default(),
    )?;

    pretty_assert_abs_diff_eq!(0.2_f64.ln(), actual.log_mass()?, epsilon = 1e-15);
    for t in [0.1, 0.01, 0.013, 1.187, 1.19] {
      pretty_assert_abs_diff_eq!(box_mass(t), (-actual.eval(t)?).exp(), epsilon = 1e-15);
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::box_narrower_than_the_support( (0.0, 1.0), 2,  0.0,     (0.0, 0.2 ), 0.0     )]
  #[case::box_wider_than_the_support(    (0.0, 1.0), 5,  LN_2,    (0.0, 2.0 ), 4.0_f64.ln())]
  #[case::shifted_box_and_support(       (0.5, 2.0), 4,  -LN_2,   (1.5, 2.25), 0.5     )]
  #[trace]
  fn test_convolution_range_function_piecewise_linear_result_keeps_the_mass_exactly(
    #[case] support: (f64, f64),
    #[case] n_points: usize,
    #[case] y: f64,
    #[case] range: (f64, f64),
    #[case] range_amplitude: f64,
  ) -> Result<(), Report> {
    let f = DistributionFunction::from_range_values(support, Array1::from_elem(n_points, y))?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?;
    let range_log_mass = (range.1 - range.0).ln() - range_amplitude;
    let function_log_mass = (support.1 - support.0).ln() - y;

    let actual = distribution_convolution(
      &DistributionNegLog::range(range, range_amplitude),
      &DistributionNegLog::Function(f),
      MaxGridPoints::default(),
    )?;

    pretty_assert_abs_diff_eq!(range_log_mass + function_log_mass, actual.log_mass()?, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_convolution_range_function_mass_converges_at_second_order() -> Result<(), Report> {
    let (lambda, range) = (1.5_f64, (0.0, 1.0));
    let errors: Vec<f64> = [32, 64, 128, 256]
      .iter()
      .map(|&cells| {
        let f = function_on((0.0, 2.0), cells + 1, |t| lambda * t)?
          .with_left_extrap(BoundaryBehavior::Hard)?
          .with_right_extrap(BoundaryBehavior::Hard)?;
        let function_mass = DistributionNegLog::Function(f.clone()).log_mass()?.exp();
        let actual = distribution_convolution(
          &DistributionNegLog::range(range, 0.0),
          &DistributionNegLog::Function(f),
          MaxGridPoints::default(),
        )?;
        Ok(actual.log_mass()?.exp() - (range.1 - range.0) * function_mass)
      })
      .try_collect::<_, _, Report>()?;
    assert_second_order(&errors);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::infinite_edge_ordinates( (0.0, 4.0), array![f64::INFINITY, 0.0, 0.0, f64::INFINITY, f64::INFINITY], vec![(0.0, 3.0)],             Array1::linspace(0.0, 3.2, 1 + 128))]
  #[case::gap_wider_than_box(      (0.0, 5.0), array![0.0, 0.0, f64::INFINITY, f64::INFINITY, 0.0, 0.0],      vec![(0.0, 2.0), (3.0, 5.0)], Array1::linspace(0.0, 5.2, 1 + 208))]
  #[trace]
  fn test_convolution_range_function_resolves_every_support_component(
    #[case] range: (f64, f64),
    #[case] y: Array1<f64>,
    #[case] support: Vec<(f64, f64)>,
    #[case] expected_t: Array1<f64>,
  ) -> Result<(), Report> {
    let f = DistributionFunction::from_range_values(range, y)?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?;

    let actual = distribution_convolution(
      &DistributionNegLog::range((0.0, 0.2), 0.0),
      &DistributionNegLog::Function(f),
      MaxGridPoints::default(),
    )?;

    let DistributionNegLog::Function(g) = &actual else {
      panic!("expected a Function result, got {actual:?}");
    };
    let overlap = |t: f64| {
      support
        .iter()
        .map(|&(s0, s1)| (t.min(s1) - (t - 0.2).max(s0)).max(0.0))
        .sum::<f64>()
    };
    pretty_assert_ulps_eq!(expected_t, g.t(), max_ulps = 4);
    pretty_assert_abs_diff_eq!(expected_t.mapv(overlap), g.y().mapv(|y| (-y).exp()), epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_convolution_range_negated_branch_ends_at_range_end() -> Result<(), Report> {
    let t = Array1::linspace(0.1, 2.1, 21);
    let branch = DistributionFunction::from_range_values((0.1, 2.1), t.mapv(|t: f64| t - 2.0 * t.ln()))?
      .with_left_extrap(BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 2.0 }))?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 }))?;
    let negated = DistributionNegLog::Function(branch).negate()?;
    let r = DistributionNegLog::range((2000.0, 2001.0), 0.0);

    let actual = distribution_convolution(&r, &negated, MaxGridPoints::default())?;

    let DistributionNegLog::Function(g) = &actual else {
      panic!("expected a Function result, got {actual:?}");
    };
    let y = g.y();
    let last = y.len() - 1;
    pretty_assert_ulps_eq!(2001.0, g.x_max(), max_ulps = 0);
    assert!(y[last].is_infinite(), "zero density at the range end, got {}", y[last]);
    assert!(
      y.slice(s![..last]).iter().all(|y| y.is_finite()),
      "density before the range end: {y:?}"
    );
    assert_eq!(BoundaryBehavior::Hard, g.right_extrap());
    let inside_last_cell = (-actual.eval(2001.0 - g.dx() / 2.0)?).exp();
    assert!(
      inside_last_cell > 0.0,
      "density {inside_last_cell} inside the last cell"
    );
    Ok(())
  }

  #[test]
  fn test_convolution_function_function_basic() -> Result<(), Report> {
    let a = DistributionNegLog::function(array![0.0, 1.0, 2.0], neg_log(&array![1.0, 2.0, 1.0]))?;
    let b = DistributionNegLog::function(array![0.0, 1.0], neg_log(&array![1.0, 2.0]))?;

    let actual = distribution_convolution(&a, &b, MaxGridPoints::default())?;

    assert_eq!(array![0.0, 1.0, 2.0, 3.0], actual.t());
    pretty_assert_abs_diff_eq!(neg_log(&array![0.25, 2.0, 2.5, 0.5]), actual.y()?, epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_convolution_single_point_functions_are_points() -> Result<(), Report> {
    let a = DistributionNegLog::function(array![2.0], array![3.0])?;
    let b = DistributionNegLog::function(array![5.0], array![4.0])?;
    let actual = distribution_convolution(&a, &b, MaxGridPoints::default())?;
    assert_eq!(DistributionNegLog::point(7.0, 7.0), actual);
    Ok(())
  }

  #[test]
  fn test_convolution_function_function_empty() -> Result<(), Report> {
    let a = DistributionNegLog::function(array![], array![])?;
    let b = DistributionNegLog::function(array![1.0, 2.0], array![0.0, 0.0])?;
    let actual = distribution_convolution(&a, &b, MaxGridPoints::default())?;
    assert_eq!(DistributionNegLog::empty(), actual);
    Ok(())
  }

  #[test]
  fn test_convolution_disjoint_grids_is_not_empty() -> Result<(), Report> {
    let a = DistributionNegLog::function(array![0.0, 1.0], array![0.0, 0.5])?;
    let b = DistributionNegLog::function(array![10.0, 11.0], array![0.0, 0.5])?;
    let actual = distribution_convolution(&a, &b, MaxGridPoints::default())?;
    assert!(!matches!(actual, DistributionNegLog::Empty));
    Ok(())
  }

  #[test]
  fn test_convolution_function_function_different_spacing() -> Result<(), Report> {
    let a = DistributionNegLog::function(array![0.0, 0.5, 1.0], neg_log(&array![1.0, 2.0, 1.0]))?;
    let b = DistributionNegLog::function(array![0.0, 1.0, 2.0], array![0.0, 0.0, 0.0])?;

    let actual = distribution_convolution(&a, &b, MaxGridPoints::default())?;

    assert_eq!(array![0.0, 1.0, 2.0, 3.0], actual.t());
    pretty_assert_abs_diff_eq!(neg_log(&array![0.125, 1.5, 1.5, 0.125]), actual.y()?, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_convolution_point_function_zero_width() -> Result<(), Report> {
    let a = DistributionNegLog::point(5.0, 1.0);
    let b = DistributionNegLog::function(array![1.0, 2.0], array![2.0, 3.0])?;

    let actual = distribution_convolution(&a, &b, MaxGridPoints::default())?;

    assert_eq!(array![6.0, 7.0], actual.t());
    pretty_assert_abs_diff_eq!(array![3.0, 4.0], actual.y()?, epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_backward_pass_temporal_direction() -> Result<(), Report> {
    let child_time_dist = DistributionNegLog::point(2013.0, 0.0);
    let negated_branch = DistributionNegLog::point(2.5, 0.0).negate()?;
    let actual = distribution_convolution(&child_time_dist, &negated_branch, MaxGridPoints::default())?;
    assert_eq!(DistributionNegLog::point(2010.5, 0.0), actual);
    Ok(())
  }

  #[test]
  fn test_forward_pass_temporal_direction() -> Result<(), Report> {
    let parent_time_dist = DistributionNegLog::point(2010.0, 0.0);
    let branch_length_dist = DistributionNegLog::point(1.5, 0.0);
    let actual = distribution_convolution(&parent_time_dist, &branch_length_dist, MaxGridPoints::default())?;
    assert_eq!(DistributionNegLog::point(2011.5, 0.0), actual);
    Ok(())
  }

  #[test]
  fn test_convolution_with_uncertainty() -> Result<(), Report> {
    let parent = DistributionNegLog::function(array![2010.0, 2010.5, 2011.0], neg_log(&array![0.2, 0.6, 0.2]))?;
    let branch = DistributionNegLog::function(array![1.0, 1.5, 2.0], neg_log(&array![0.3, 0.4, 0.3]))?;

    let actual = distribution_convolution(&parent, &branch, MaxGridPoints::default())?;

    assert_eq!(array![2011.0, 2011.5, 2012.0, 2012.5, 2013.0], actual.t());
    pretty_assert_abs_diff_eq!(
      neg_log(&array![0.0075, 0.065, 0.15, 0.065, 0.0075]),
      actual.y()?,
      epsilon = 1e-12
    );
    Ok(())
  }

  #[test]
  fn test_convolution_convolve_small_dx_function_function() -> Result<(), Report> {
    let dx = 1e-7;
    let num_points = 500;
    let values = Array1::zeros(num_points);
    let dist_a = DistributionNegLog::Function(DistributionFunction::from_start_dx_values(0.0, dx, values.clone())?);
    let dist_b = DistributionNegLog::Function(DistributionFunction::from_start_dx_values(0.0, dx, values)?);

    let actual = distribution_convolution(&dist_a, &dist_b, MaxGridPoints::default())?;

    let DistributionNegLog::Function(result_fn) = actual else {
      panic!("expected a Function result, got {actual:?}");
    };
    assert_eq!(2 * num_points - 1, result_fn.len());
    pretty_assert_abs_diff_eq!(0.0, result_fn.x_min(), epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(dx, result_fn.dx(), epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_convolution_convolve_small_dx_range_function() -> Result<(), Report> {
    let dx = 1e-7;
    let n = 500;
    let func = DistributionNegLog::Function(
      DistributionFunction::from_start_dx_values(0.0, dx, Array1::zeros(n))?
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Hard)?,
    );
    let range = DistributionNegLog::range((0.0, 2e-7), 0.0);

    let actual = distribution_convolution(&range, &func, MaxGridPoints::default())?;

    let DistributionNegLog::Function(f) = actual else {
      panic!("expected a Function result, got {actual:?}");
    };
    let rise_cell = 2e-7 / 8.0;
    assert_eq!(1 + (n + 1) * 4, f.len());
    pretty_assert_abs_diff_eq!(rise_cell, f.dx(), epsilon = 1e-22);
    Ok(())
  }

  #[test]
  fn test_coarsen_convolution_narrow_range_keeps_fine_grid() -> Result<(), Report> {
    let fine = DistributionFunction::from_start_dx_values(0.0, 0.01, array![1.0, 0.0, 1.0])?;
    assert_error!(
      Grid::from_range_dx(0.0, 0.02, 1.0, MaxGridPoints::default()),
      "Grid must have at least 2 points, got 1"
    );

    let actual = coarsen_convolution(fine.clone(), 1.0, MaxGridPoints::default())?;

    assert_eq!(DistributionNegLog::Function(fine), actual);
    Ok(())
  }

  #[test]
  fn test_coarsen_convolution_half_cell_range_resamples_to_two_points() -> Result<(), Report> {
    let fine = DistributionFunction::from_start_dx_values(0.0, 0.25, array![1.0, 0.0, 1.0])?;
    let expected_points = Grid::from_range_dx(0.0, 0.5, 1.0, MaxGridPoints::default())?.n_points();
    assert_eq!(2, expected_points);

    let DistributionNegLog::Function(actual) = coarsen_convolution(fine, 1.0, MaxGridPoints::default())? else {
      panic!("expected a Function result");
    };

    assert_eq!(expected_points, actual.len());
    pretty_assert_abs_diff_eq!(1.0, actual.dx(), epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_coarsen_convolution_wide_range_resamples_to_coarse_spacing() -> Result<(), Report> {
    let fine = DistributionFunction::from_start_dx_values(0.0, 0.5, Array1::zeros(9))?;
    let expected_points = Grid::from_range_dx(0.0, 4.0, 1.0, MaxGridPoints::default())?.n_points();
    assert_eq!(5, expected_points);

    let DistributionNegLog::Function(actual) = coarsen_convolution(fine, 1.0, MaxGridPoints::default())? else {
      panic!("expected a Function result");
    };

    assert_eq!(expected_points, actual.len());
    pretty_assert_abs_diff_eq!(1.0, actual.dx(), epsilon = 1e-15);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::distribution_core::formula::DistributionFormula;
    use crate::policy::NegLog;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum DistributionVariant {
      Empty,
      Point,
      Range,
      Function,
      Formula,
    }

    pub(super) fn distribution(variant: DistributionVariant) -> DistributionNegLog {
      match variant {
        DistributionVariant::Empty => DistributionNegLog::empty(),
        DistributionVariant::Point => DistributionNegLog::point(0.0, 0.0),
        DistributionVariant::Range => DistributionNegLog::range((0.0, 1.0), 0.0),
        DistributionVariant::Function => {
          DistributionNegLog::function(array![0.0, 1.0, 2.0], array![1.0, 0.0, 1.0]).unwrap()
        },
        DistributionVariant::Formula => DistributionNegLog::Formula(DistributionFormula::new(|_| Ok(0.0), 0.0, 1.0)),
      }
    }

    pub(super) fn function_on(
      range: (f64, f64),
      n_points: usize,
      y: impl Fn(f64) -> f64,
    ) -> Result<DistributionFunction<f64, NegLog>, Report> {
      DistributionFunction::from_range_values(range, Array1::linspace(range.0, range.1, n_points).mapv(y))
    }

    pub(super) fn assert_grid_point(actual: &DistributionNegLog, t: f64) {
      let DistributionNegLog::Function(g) = actual else {
        panic!("expected a Function result, got {actual:?}");
      };
      assert!(
        g.t().iter().any(|&x| ulps_eq!(x, t, epsilon = 0.0, max_ulps = 4)),
        "{t} is a grid point of the result"
      );
    }

    pub(super) fn neg_log(plain: &Array1<f64>) -> Array1<f64> {
      plain.mapv(|p| -p.ln())
    }

    pub(super) fn exponential_box_density(lambda: f64, (a, b): (f64, f64), t: f64) -> f64 {
      if t <= a {
        0.0
      } else if t <= b {
        -(-lambda * (t - a)).exp_m1()
      } else {
        (-lambda * (t - b)).exp() * -(-lambda * (b - a)).exp_m1()
      }
    }

    pub(super) struct PiecewiseLinear {
      x: Array1<f64>,
      p: Array1<f64>,
      right_rate: f64,
    }

    impl PiecewiseLinear {
      pub(super) fn of(f: &DistributionFunction<f64, NegLog>, right_rate: f64) -> Self {
        Self {
          x: f.t(),
          p: f.y().mapv(|y| (-y).exp()),
          right_rate,
        }
      }

      pub(super) fn box_mass(&self, (start, end): (f64, f64), time: f64) -> f64 {
        let (from, to) = (time - end, time - start);
        let last = self.x.len() - 1;
        let cells: f64 = (0..last)
          .map(|cell| {
            let (x0, x1) = (self.x[cell], self.x[cell + 1]);
            let (lower, upper) = (from.max(x0), to.min(x1));
            if lower >= upper {
              return 0.0;
            }
            let density = |at: f64| self.p[cell] + (self.p[cell + 1] - self.p[cell]) * (at - x0) / (x1 - x0);
            (upper - lower) * f64::midpoint(density(lower), density(upper))
          })
          .sum();
        let edge = self.x[last];
        let (near, far) = (from.max(edge) - edge, to.max(edge) - edge);
        let tail = if self.right_rate.is_finite() {
          self.p[last] * ((-self.right_rate * near).exp() - (-self.right_rate * far).exp()) / self.right_rate
        } else {
          0.0
        };
        cells + tail
      }
    }
  }
}
