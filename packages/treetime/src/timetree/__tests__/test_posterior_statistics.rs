#[cfg(test)]
mod tests {
  use crate::timetree::max_posterior_region::{PosteriorSurvival, posterior_median, posterior_survival};
  use eyre::Report;
  use helpers::{flat_unit_function, grid};
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use statrs::distribution::{ContinuousCDF, Gamma};
  use treetime_distribution::{BoundaryBehavior, Distribution, DistributionFunction, NegLog};
  use treetime_grid::{HardApproachLaw, SoftTailLaw};
  use treetime_utils::pretty_assert_abs_diff_eq;

  const TRIM: f64 = 1e-9;
  const MAX_ERROR: f64 = 1e-6;

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
  fn test_posterior_survival_of_point_and_range() -> Result<(), Report> {
    let point = PosteriorSurvival {
      t: vec![2015.5],
      survival: vec![0.0],
    };
    let range = PosteriorSurvival {
      t: vec![2015.0, 2017.0],
      survival: vec![1.0, 0.0],
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
    use crate::timetree::max_posterior_region::PosteriorSurvival;
    use eyre::Report;
    use itertools::Itertools;
    use ndarray::{Array1, array};
    use treetime_distribution::{BoundaryBehavior, Distribution, DistributionFunction, NegLog};
    use treetime_utils::cast::ToF64Exact;
    use treetime_utils::pretty_assert_abs_diff_eq;

    pub(super) fn assert_interpolates_within_max_error(survival: &PosteriorSurvival, exact: impl Fn(f64) -> f64) {
      let PosteriorSurvival { t, survival } = survival;
      assert!(t.iter().tuple_windows().all(|(a, b)| a < b));
      for (i, (t0, t1)) in t.iter().tuple_windows().enumerate() {
        for k in 0..=16 {
          let x = t0 + (t1 - t0) * (k.to_f64_exact() / 16.0);
          let interpolated = survival[i] + (x - t0) / (t1 - t0) * (survival[i + 1] - survival[i]);
          pretty_assert_abs_diff_eq!(exact(x), interpolated, epsilon = 1e-6);
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
