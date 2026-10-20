#[cfg(test)]
mod tests {
  use crate::coalescent::mean_field::NodeSurvival;
  use crate::timetree::max_posterior_region::{PosteriorSurvival, posterior_median, posterior_survival};
  use eyre::Report;
  use helpers::{flat_unit_function, grid};
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use statrs::distribution::{ContinuousCDF, Gamma};
  use treetime_distribution::{BoundaryBehavior, Distribution, DistributionFunction, MassProfile, NegLog, TailMass};
  use treetime_grid::{HardApproachLaw, SoftTailLaw};
  use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq};

  const TRIM: f64 = 1e-9;
  const MAX_ERROR: f64 = 1e-6;
  const CALENDAR_TIME: f64 = 2000.0;
  const STEP: f64 = 1.0 / 1_048_576.0;

  #[test]
  fn test_posterior_median_of_a_skewed_gamma() -> Result<(), Report> {
    let (shape, rate) = (3.0, 2.0);
    let t = grid(1e-4, 12.0, 120_000);
    let neglog = t.mapv(|x| 5.0 - (shape - 1.0) * x.ln() + rate * x);

    let median = posterior_median(&Distribution::<NegLog>::function(t, neglog)?)?;

    pretty_assert_abs_diff_eq!(
      Gamma::new(shape, rate)?.inverse_cdf(0.5),
      median.unwrap(),
      epsilon = 1e-7
    );
    Ok(())
  }

  #[test]
  fn test_posterior_median_lies_between_two_equal_peaks() -> Result<(), Report> {
    let t = grid(2010.0, 2018.0, 80_001);
    let neglog = t.mapv(|x| {
      let a = (-(x - 2013.0).powi(2) / 0.02).exp();
      let b = (-(x - 2015.0).powi(2) / 0.02).exp();
      -(a + b).ln()
    });

    let median = posterior_median(&Distribution::<NegLog>::function(t, neglog)?)?.unwrap();

    assert!(median > 2013.5 && median < 2014.5, "median {median}");
    Ok(())
  }

  #[test]
  fn test_posterior_median_of_point_range_and_empty() -> Result<(), Report> {
    assert_eq!(
      Some(2015.5),
      posterior_median(&Distribution::<NegLog>::point(2015.5, 0.0))?
    );
    assert_eq!(
      Some(2016.0),
      posterior_median(&Distribution::<NegLog>::range((2015.0, 2017.0), 3.0))?
    );
    assert_eq!(None, posterior_median(&Distribution::<NegLog>::Empty)?);
    Ok(())
  }

  #[test]
  fn test_posterior_median_counts_mass_of_the_earlier_tail() -> Result<(), Report> {
    let unit_box_with_tail_mass_one = flat_unit_function(
      (-1.0, 0.0),
      BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 }),
      BoundaryBehavior::Hard,
    )?;

    let median = posterior_median(&unit_box_with_tail_mass_one)?;

    assert_eq!(Some(-1.0), median);
    Ok(())
  }

  #[test]
  fn test_posterior_survival_continues_into_the_earlier_tail_up_to_the_trim() -> Result<(), Report> {
    let unit_box_with_tail_mass_one = flat_unit_function(
      (-1.0, 0.0),
      BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 }),
      BoundaryBehavior::Hard,
    )?;

    let survival = posterior_survival(&unit_box_with_tail_mass_one, TRIM, MAX_ERROR)?.unwrap();

    let exact = |t: f64| {
      if t < -1.0 {
        1.0 - (t + 1.0).exp() / 2.0
      } else {
        -t / 2.0
      }
    };
    let last = survival.t.len() - 1;
    pretty_assert_abs_diff_eq!(-1.0 + (2.0 * TRIM).ln(), survival.t[0], epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(1.0 - TRIM, survival.survival[0], epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(-2.0 * TRIM, survival.t[last], epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(TRIM, survival.survival[last], epsilon = 1e-15);
    helpers::assert_interpolates_within_max_error(&survival, exact);
    Ok(())
  }

  #[test]
  fn test_posterior_survival_follows_a_steep_later_tail() -> Result<(), Report> {
    let rate = 100.0;
    let distribution = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values((0.0, 1.0), Array1::zeros(11))?
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: rate }))?,
    );

    let survival = posterior_survival(&distribution, TRIM, MAX_ERROR)?.unwrap();

    let total = 1.0 + 1.0 / rate;
    let exact = |t: f64| {
      if t <= 1.0 {
        (1.0 - t + 1.0 / rate) / total
      } else {
        (-rate * (t - 1.0)).exp() / (rate * total)
      }
    };
    let tail_end = 1.0 + (1.0 / (rate * total * TRIM)).ln() / rate;
    let last = survival.t.len() - 1;
    pretty_assert_abs_diff_eq!(TRIM * total, survival.t[0], epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(tail_end, survival.t[last], epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(TRIM, survival.survival[last], epsilon = 1e-20);
    helpers::assert_interpolates_within_max_error(&survival, exact);
    Ok(())
  }

  #[test]
  fn test_posterior_survival_ends_inside_a_hard_approach_tail() -> Result<(), Report> {
    let distribution = flat_unit_function(
      (0.0, 1.0),
      BoundaryBehavior::Hard,
      BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 2.0, b: 1.0 }),
    )?;

    let survival = posterior_survival(&distribution, TRIM, MAX_ERROR)?.unwrap();

    let total = 1.5;
    let exact = |t: f64| {
      if t <= 1.0 {
        (1.5 - t) / total
      } else {
        (2.0 - t).powi(2) / (2.0 * total)
      }
    };
    let last = survival.t.len() - 1;
    pretty_assert_abs_diff_eq!(2.0 - (2.0 * total * TRIM).sqrt(), survival.t[last], epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(TRIM, survival.survival[last], epsilon = 1e-20);
    assert!(survival.t[last] < 2.0);
    helpers::assert_interpolates_within_max_error(&survival, exact);
    Ok(())
  }

  #[test]
  fn test_posterior_survival_is_exact_inside_regular_cells_with_knots_only_at_grid_points() -> Result<(), Report> {
    let density = [1.0, 3.0, 2.0, 4.0, 1.0];
    let distribution = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values(
        (0.0, 4.0),
        Array1::from_iter(density.map(|p: f64| -p.ln())),
      )?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?,
    );

    let survival = posterior_survival(&distribution, TRIM, MAX_ERROR)?.unwrap();

    let cell_mass_above = |cell: usize, t: f64| {
      let at = density[cell] + (density[cell + 1] - density[cell]) * (t - cell.to_f64_exact());
      (cell.to_f64_exact() + 1.0 - t) * f64::midpoint(at, density[cell + 1])
    };
    let exact = |t: f64| {
      let cell = f64_to_usize_saturating(t.floor()).min(3);
      let later: f64 = (cell + 1..4)
        .map(|later| cell_mass_above(later, later.to_f64_exact()))
        .sum();
      (cell_mass_above(cell, t) + later) / 10.0
    };
    assert_eq!(vec![1.0, 2.0, 3.0], survival.t[1..4].to_vec());
    assert_eq!(5, survival.t.len());
    helpers::assert_interpolates_to_rounding(&survival, exact);
    Ok(())
  }

  #[test]
  fn test_posterior_survival_is_exact_across_flat_hard_cells_around_an_interior_gap() -> Result<(), Report> {
    let distribution = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values(
        (0.0, 5.0),
        array![0.0, 0.0, f64::INFINITY, f64::INFINITY, 0.0, 0.0],
      )?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?,
    );

    let survival = posterior_survival(&distribution, TRIM, MAX_ERROR)?.unwrap();

    let exact = |t: f64| ((2.0 - t).max(0.0) + (5.0 - t.max(3.0))) / 4.0;
    assert_eq!(vec![1.0, 2.0, 2.0, 3.0, 3.0, 4.0], survival.t[1..7].to_vec());
    helpers::assert_interpolates_to_rounding(&survival, exact);
    Ok(())
  }

  #[test]
  fn test_posterior_survival_follows_the_power_law_of_a_hard_cell_end() -> Result<(), Report> {
    let (anchor_y, b) = (0.4, 2.5);
    let distribution = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values(
        (0.0, 3.0),
        array![anchor_y, anchor_y, anchor_y, f64::INFINITY],
      )?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 3.0, b }))?,
    );

    let survival = posterior_survival(&distribution, TRIM, MAX_ERROR)?.unwrap();

    let power = b + 1.0;
    let total = 2.0 + 1.0 / power;
    let exact = |t: f64| {
      if t <= 2.0 {
        (2.0 - t + 1.0 / power) / total
      } else {
        (3.0 - t).powf(power) / power / total
      }
    };
    let last = survival.t.len() - 1;
    pretty_assert_abs_diff_eq!(
      3.0 - (power * total * TRIM).powf(power.recip()),
      survival.t[last],
      epsilon = 1e-14
    );
    helpers::assert_interpolates_within_max_error(&survival, exact);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::below_linear(   0.5)]
  #[case::between_2_and_3(2.5)]
  #[trace]
  fn test_posterior_survival_follows_the_power_law_at_hard_points_inside_the_window(#[case] b: f64) -> Result<(), Report> {
    let edge = -b * 2.0_f64.ln();
    let distribution = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values(
        (0.0, 5.0),
        array![edge, 0.0, f64::INFINITY, f64::INFINITY, 0.0, edge],
      )?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?,
    );

    let survival = posterior_survival(&distribution, TRIM, MAX_ERROR)?.unwrap();

    let (p_edge, power) = (2.0_f64.powf(b), b + 1.0);
    let trapezoid = f64::midpoint(p_edge, 1.0);
    let total = 2.0 * (trapezoid + 1.0 / power);
    let mass_above = |t: f64| {
      if t <= 1.0 {
        let u = 1.0 - t;
        u * f64::midpoint(1.0 + (p_edge - 1.0) * u, 1.0) + 1.0 / power + total / 2.0
      } else if t <= 2.0 {
        (2.0 - t).powf(power) / power + total / 2.0
      } else if t <= 3.0 {
        total / 2.0
      } else if t <= 4.0 {
        trapezoid + (1.0 - (t - 3.0).powf(power)) / power
      } else {
        let u = 5.0 - t;
        u * f64::midpoint(p_edge + (1.0 - p_edge) * u, p_edge)
      }
    };
    let near_hard_points = (0..=40).flat_map(|k| {
      let distance = 10.0_f64.powf(-k.to_f64_exact() / 4.0);
      [2.0 - distance, 2.0 + distance, 3.0 - distance, 3.0 + distance]
    });
    let interpolant = NodeSurvival::from(survival.clone());
    for t in near_hard_points.filter(|&t| t > 1.0 && t < 4.0) {
      pretty_assert_abs_diff_eq!(mass_above(t) / total, interpolant.at(t), epsilon = 1e-6);
    }
    helpers::assert_interpolates_within_max_error(&survival, |t| mass_above(t) / total);
    Ok(())
  }

  #[test]
  fn test_posterior_survival_follows_a_long_soft_tail_of_a_narrow_posterior() -> Result<(), Report> {
    let distribution = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values((0.0, 0.01), Array1::zeros(11))?
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 }))?,
    );

    let survival = posterior_survival(&distribution, TRIM, MAX_ERROR)?.unwrap();

    let total = 1.01;
    let exact = |t: f64| {
      if t <= 0.01 {
        (1.01 - t) / total
      } else {
        (-(t - 0.01)).exp() / total
      }
    };
    let last = survival.t.len() - 1;
    pretty_assert_abs_diff_eq!(0.01 + (1.0 / (total * TRIM)).ln(), survival.t[last], epsilon = 1e-13);
    helpers::assert_interpolates_within_max_error(&survival, exact);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::earlier_hard_approach( CALENDAR_TIME + STEP, array![0.0, 0.0],           BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: CALENDAR_TIME, b: 0.0 }), BoundaryBehavior::Hard)]
  #[case::later_hard_approach(   CALENDAR_TIME,        array![0.0, 0.0],           BoundaryBehavior::Hard, BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: CALENDAR_TIME + 2.0 * STEP, b: 0.0 }))]
  #[case::interior_hard_point(   CALENDAR_TIME,        array![0.0, f64::INFINITY, 0.0], BoundaryBehavior::Hard, BoundaryBehavior::Hard)]
  #[trace]
  fn test_node_survival_of_a_uniform_density_that_reaches_a_flat_hard_point_is_exact(
    #[case] start: f64,
    #[case] y: Array1<f64>,
    #[case] left: BoundaryBehavior,
    #[case] right: BoundaryBehavior,
  ) -> Result<(), Report> {
    let distribution = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_start_dx_values(start, STEP, y)?
        .with_left_extrap(left)?
        .with_right_extrap(right)?,
    );

    let survival = NodeSurvival::of_posterior(&distribution)?.unwrap();

    let end = CALENDAR_TIME + 2.0 * STEP;
    for k in 0..=64 {
      let t = CALENDAR_TIME + 2.0 * STEP * (k.to_f64_exact() / 64.0);
      pretty_assert_abs_diff_eq!((end - t) / (2.0 * STEP), survival.at(t), epsilon = 1e-8);
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::flat_cell(       array![0.0, 0.0]               )]
  #[case::trapezoid_cells( array![0.0, 0.5, 0.2, 1.0, 0.3])]
  #[trace]
  fn test_node_survival_follows_the_mass_model_on_cells_a_few_hundred_coordinate_ulps_wide(
    #[case] y: Array1<f64>,
  ) -> Result<(), Report> {
    let function = DistributionFunction::<f64, NegLog>::from_start_dx_values(CALENDAR_TIME, 1e-10, y)?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?;
    let profile = MassProfile::new(&function)?;
    let fractions = profile.fractions(TailMass::Required).unwrap();

    let survival = NodeSurvival::of_posterior(&Distribution::Function(function.clone()))?.unwrap();

    let samples = 64 * (function.y().len() - 1);
    for k in 0..=samples {
      let t = function.x_min() + (function.x_max() - function.x_min()) * (k.to_f64_exact() / samples.to_f64_exact());
      pretty_assert_abs_diff_eq!(fractions.survival(t), survival.at(t), epsilon = 1e-8);
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::right((-1.0, 0.0), array![0.0, -(1e-300_f64).ln()], BoundaryBehavior::Hard, BoundaryBehavior::Linear(SoftTailLaw { slope: 1e-310 }), "Posterior on [-1, 0] has a tail of slope 1e-310 on the Right side, which holds mass 1e-9 beyond the largest finite time")]
  #[case::left( ( 0.0, 1.0), array![-(1e-300_f64).ln(), 0.0], BoundaryBehavior::Linear(SoftTailLaw { slope: -1e-310 }), BoundaryBehavior::Hard, "Posterior on [0, 1] has a tail of slope -1e-310 on the Left side, which holds mass 1e-9 beyond the largest finite time")]
  #[trace]
  fn test_posterior_survival_of_a_tail_too_slow_for_finite_times_is_error(
    #[case] range: (f64, f64),
    #[case] y: Array1<f64>,
    #[case] left: BoundaryBehavior,
    #[case] right: BoundaryBehavior,
    #[case] message: &str,
  ) -> Result<(), Report> {
    let distribution = Distribution::Function(
      DistributionFunction::<f64, NegLog>::from_range_values(range, y)?
        .with_left_extrap(left)?
        .with_right_extrap(right)?,
    );
    assert_error!(posterior_survival(&distribution, TRIM, MAX_ERROR), message);
    Ok(())
  }

  #[test]
  fn test_posterior_survival_of_point_and_range() -> Result<(), Report> {
    let point = PosteriorSurvival {
      t: vec![2015.5],
      survival: vec![0.0],
      slope: vec![0.0],
    };
    let range = PosteriorSurvival {
      t: vec![2015.0, 2017.0],
      survival: vec![1.0, 0.0],
      slope: vec![-0.5, -0.5],
    };
    assert_eq!(
      Some(point),
      posterior_survival(&Distribution::<NegLog>::point(2015.5, 0.0), TRIM, MAX_ERROR)?
    );
    assert_eq!(
      Some(range),
      posterior_survival(&Distribution::<NegLog>::range((2015.0, 2017.0), 3.0), TRIM, MAX_ERROR)?
    );
    assert_eq!(
      None,
      posterior_survival(&Distribution::<NegLog>::Empty, TRIM, MAX_ERROR)?
    );
    Ok(())
  }

  mod helpers {
    use crate::coalescent::mean_field::NodeSurvival;
    use crate::timetree::max_posterior_region::PosteriorSurvival;
    use eyre::Report;
    use itertools::Itertools;
    use ndarray::{Array1, array};
    use treetime_distribution::{BoundaryBehavior, Distribution, DistributionFunction, NegLog};
    use treetime_utils::cast::ToF64Exact;
    use treetime_utils::pretty_assert_abs_diff_eq;

    pub(super) fn assert_interpolates_within_max_error(survival: &PosteriorSurvival, exact: impl Fn(f64) -> f64) {
      for_each_interpolated(survival, |x, interpolated| {
        pretty_assert_abs_diff_eq!(exact(x), interpolated, epsilon = 1e-6);
      });
    }

    pub(super) fn assert_interpolates_to_rounding(survival: &PosteriorSurvival, exact: impl Fn(f64) -> f64) {
      for_each_interpolated(survival, |x, interpolated| {
        pretty_assert_abs_diff_eq!(exact(x), interpolated, epsilon = 1e-15);
      });
    }

    fn for_each_interpolated(survival: &PosteriorSurvival, mut check: impl FnMut(f64, f64)) {
      let t = &survival.t;
      assert!(t.iter().tuple_windows().all(|(a, b)| a <= b));
      let interpolant = NodeSurvival::from(survival.clone());
      for (t0, t1) in t.iter().tuple_windows() {
        for k in 0..16 {
          let x = t0 + (t1 - t0) * (k.to_f64_exact() / 16.0);
          check(x, interpolant.at(x));
        }
      }
    }

    pub(super) fn grid(start: f64, end: f64, n_points: usize) -> Array1<f64> {
      Array1::linspace(start, end, n_points)
    }

    pub(super) fn flat_unit_function(
      range: (f64, f64),
      left: BoundaryBehavior,
      right: BoundaryBehavior,
    ) -> Result<Distribution<NegLog>, Report> {
      let function = DistributionFunction::from_range_values(range, array![0.0, 0.0])?
        .with_left_extrap(left)?
        .with_right_extrap(right)?;
      Ok(Distribution::Function(function))
    }
  }
}
