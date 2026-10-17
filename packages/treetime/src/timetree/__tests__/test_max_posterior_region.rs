#[cfg(test)]
mod tests {
  use crate::timetree::max_posterior_region::{
    PeakLocation, PosteriorRegion, grid_max_posterior_region, max_posterior_region, posterior_median,
  };
  use eyre::Report;
  use helpers::{grid, interpolate};
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use statrs::distribution::{Continuous, ContinuousCDF, Gamma, Normal};
  use treetime_distribution::{Distribution, NegLog};
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq, prop_assert_abs_diff_eq};

  #[test]
  fn test_max_posterior_region_gaussian_is_symmetric_hpd() -> Result<(), Report> {
    let (mean, std_dev) = (2015.0, 0.5);
    let t = grid(mean - 8.0 * std_dev, mean + 8.0 * std_dev, 160_001);
    let neglog = t.mapv(|x| 3.7 + (x - mean).powi(2) / (2.0 * std_dev * std_dev));

    let region = grid_max_posterior_region(&t, &neglog, 0.9)?;

    let normal = Normal::new(mean, std_dev)?;
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_abs_diff_eq!(normal.inverse_cdf(0.05), region.lower, epsilon = 1e-8);
    pretty_assert_abs_diff_eq!(normal.inverse_cdf(0.95), region.upper, epsilon = 1e-8);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_skewed_gamma_has_equal_density_bounds() -> Result<(), Report> {
    let (shape, rate) = (3.0, 2.0);
    let t = grid(1e-4, 12.0, 120_000);
    let neglog = t.mapv(|x| -(shape - 1.0) * x.ln() + rate * x);

    let region = grid_max_posterior_region(&t, &neglog, 0.9)?;

    let gamma = Gamma::new(shape, rate)?;
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_abs_diff_eq!(0.9, gamma.cdf(region.upper) - gamma.cdf(region.lower), epsilon = 1e-8);
    pretty_assert_abs_diff_eq!(gamma.pdf(region.lower), gamma.pdf(region.upper), epsilon = 1e-8);
    Ok(())
  }

  #[test]
  fn test_posterior_median_of_a_skewed_gamma() -> Result<(), Report> {
    let (shape, rate) = (3.0, 2.0);
    let t = grid(1e-4, 12.0, 120_000);
    let neglog = t.mapv(|x| 5.0 - (shape - 1.0) * x.ln() + rate * x);

    let median = posterior_median(&Distribution::<NegLog>::function(t, neglog)?)?;

    pretty_assert_abs_diff_eq!(
      Gamma::new(shape, rate)?.inverse_cdf(0.5),
      median.unwrap(),
      epsilon = 1e-6
    );
    Ok(())
  }

  #[test]
  fn test_posterior_median_lies_between_two_equal_peaks() -> Result<(), Report> {
    // two equal, narrow peaks at 2013 and 2015: the peak is either, the median is between them
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
  fn test_posterior_median_of_point_and_empty() -> Result<(), Report> {
    assert_eq!(
      Some(2015.5),
      posterior_median(&Distribution::<NegLog>::point(2015.5, 0.0))?
    );
    assert_eq!(None, posterior_median(&Distribution::<NegLog>::Empty)?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::level_below_the_second_mode_keeps_main_mode(0.6,  false)]
  #[case::level_between_second_mode_and_dip(           0.8,  false)]
  #[case::level_above_the_dip_adds_second_mode(       0.99, true)]
  #[trace]
  fn test_max_posterior_region_bimodal_region_is_level_set_connected_to_peak(
    #[case] fraction: f64,
    #[case] contains_second_mode: bool,
  ) -> Result<(), Report> {
    let main = Normal::new(5.0, 0.5)?;
    let side_mode = 3.2;
    let side = Normal::new(side_mode, 0.3)?;
    let (main_weight, side_weight) = (0.85, 0.15);
    let t = grid(-2.0, 10.0, 120_001);
    let neglog = t.mapv(|x| -(main_weight * main.pdf(x) + side_weight * side.pdf(x)).ln());

    let region = grid_max_posterior_region(&t, &neglog, fraction)?;

    let mixture_cdf = |x: f64| main_weight * main.cdf(x) + side_weight * side.cdf(x);
    let level = interpolate(&t, &neglog, region.upper);
    let highest_inside = t
      .iter()
      .zip(&neglog)
      .filter(|(x, _)| **x > region.lower && **x < region.upper)
      .map(|(_, y)| *y)
      .fold(f64::NEG_INFINITY, f64::max);
    assert_eq!(PeakLocation::Interior, region.peak);
    assert_eq!(contains_second_mode, region.lower < side_mode);
    pretty_assert_abs_diff_eq!(fraction, mixture_cdf(region.upper) - mixture_cdf(region.lower), epsilon = 1e-8);
    pretty_assert_abs_diff_eq!(level, interpolate(&t, &neglog, region.lower), epsilon = 1e-9);
    assert!(highest_inside <= level, "every point inside {region:?} lies above the density level");
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_peak_at_earliest_date_extends_later() -> Result<(), Report> {
    let t = grid(2000.0, 2010.0, 100_001);
    let neglog = t.mapv(|x| x - 2000.0);

    let region = grid_max_posterior_region(&t, &neglog, 0.9)?;

    let truncated_exponential_quantile = 2000.0 - (1.0 - 0.9 * (1.0 - (-10.0_f64).exp())).ln();
    assert_eq!(PeakLocation::Earliest, region.peak);
    pretty_assert_abs_diff_eq!(2000.0, region.lower, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(truncated_exponential_quantile, region.upper, epsilon = 1e-9);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_peak_at_latest_date_extends_earlier() -> Result<(), Report> {
    let t = grid(2000.0, 2010.0, 100_001);
    let neglog = t.mapv(|x| 2010.0 - x);

    let region = grid_max_posterior_region(&t, &neglog, 0.9)?;

    let truncated_exponential_quantile = 2010.0 + (1.0 - 0.9 * (1.0 - (-10.0_f64).exp())).ln();
    assert_eq!(PeakLocation::Latest, region.peak);
    pretty_assert_abs_diff_eq!(truncated_exponential_quantile, region.lower, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(2010.0, region.upper, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_flat_range_resolves_tie_to_latest_date() -> Result<(), Report> {
    let distribution = Distribution::<NegLog>::range((2019.0, 2021.0), 4.0);

    let region = max_posterior_region(&distribution, 0.9)?;

    let expected = PosteriorRegion {
      lower: 2019.2,
      upper: 2021.0,
      peak: PeakLocation::Latest,
    };
    assert_eq!(Some(expected), region);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_point_is_zero_width() -> Result<(), Report> {
    let distribution = Distribution::<NegLog>::point(2020.5, 1.0);

    let region = max_posterior_region(&distribution, 0.9)?;

    let expected = PosteriorRegion {
      lower: 2020.5,
      upper: 2020.5,
      peak: PeakLocation::Latest,
    };
    assert_eq!(Some(expected), region);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_zero_probability_point_has_no_region() -> Result<(), Report> {
    let distribution = Distribution::<NegLog>::point(2020.5, f64::INFINITY);
    assert_eq!(None, max_posterior_region(&distribution, 0.9)?);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_empty_has_no_region() -> Result<(), Report> {
    assert_eq!(None, max_posterior_region(&Distribution::<NegLog>::empty(), 0.9)?);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_ignores_zero_probability_tails() -> Result<(), Report> {
    let (mean, std_dev) = (2015.0, 0.5);
    let t = grid(mean - 8.0 * std_dev, mean + 8.0 * std_dev, 16001);
    let neglog = t.mapv(|x| (x - mean).powi(2) / (2.0 * std_dev * std_dev));
    let neglog_hard_tails = neglog.mapv(|y| if y > 30.0 { f64::INFINITY } else { y });

    let region = grid_max_posterior_region(&t, &neglog, 0.9)?;
    let region_hard_tails = grid_max_posterior_region(&t, &neglog_hard_tails, 0.9)?;

    assert_eq!(PeakLocation::Interior, region_hard_tails.peak);
    pretty_assert_abs_diff_eq!(region.lower, region_hard_tails.lower, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(region.upper, region_hard_tails.upper, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_max_posterior_region_rejects_nan_ordinates() {
    let result = grid_max_posterior_region(&array![0.0, 1.0, 2.0], &array![0.0, f64::NAN, 1.0], 0.9);
    assert_error!(result, "Posterior on [0, 2] has NaN negative log-likelihood values");
  }

  #[test]
  fn test_max_posterior_region_rejects_unsorted_grid() {
    let result = grid_max_posterior_region(&array![0.0, 2.0, 1.0], &array![1.0, 0.0, 1.0], 0.9);
    assert_error!(result, "Time grid of a posterior must be strictly increasing");
  }

  #[test]
  fn test_max_posterior_region_rejects_grid_without_mass() {
    let result = grid_max_posterior_region(&array![0.0, 1.0], &array![f64::INFINITY, f64::INFINITY], 0.9);
    assert_error!(result, "Posterior on [0, 1] has no probability mass");
  }

  proptest! {
    #[test]
    fn test_prop_max_posterior_region_interior_bounds_have_equal_density(
      center in 0.4_f64..0.6,
      width_lower in 0.01_f64..0.1,
      width_upper in 0.01_f64..0.1,
      power in 1.0_f64..3.0,
    ) {
      let t = grid(0.0, 1.0, 2001);
      let neglog = t.mapv(|x| {
        let width = if x < center { width_lower } else { width_upper };
        ((x - center).abs() / width).powf(power)
      });

      let region = grid_max_posterior_region(&t, &neglog, 0.9).unwrap();

      prop_assume!(region.peak == PeakLocation::Interior);
      prop_assert_abs_diff_eq!(interpolate(&t, &neglog, region.lower), interpolate(&t, &neglog, region.upper), epsilon = 1e-9);
    }

    #[test]
    fn test_prop_max_posterior_region_affine_time_change_preserves_region(
      center in 0.2_f64..0.8,
      width in 0.02_f64..0.3,
      scale in 0.1_f64..10.0,
      shift in -3000.0_f64..3000.0,
    ) {
      let t = grid(0.0, 1.0, 1001);
      let neglog = t.mapv(|x| ((x - center) / width).powi(2) + (x - center));

      let region = grid_max_posterior_region(&t, &neglog, 0.9).unwrap();
      let moved = grid_max_posterior_region(&t.mapv(|x| scale * x + shift), &neglog, 0.9).unwrap();

      prop_assert_eq!(region.peak, moved.peak);
      prop_assert_abs_diff_eq!(scale * region.lower + shift, moved.lower, epsilon = 1e-9);
      prop_assert_abs_diff_eq!(scale * region.upper + shift, moved.upper, epsilon = 1e-9);
    }
  }

  mod helpers {
    use ndarray::Array1;

    pub(super) fn grid(start: f64, end: f64, n_points: usize) -> Array1<f64> {
      Array1::linspace(start, end, n_points)
    }

    pub(super) fn interpolate(t: &Array1<f64>, y: &Array1<f64>, x: f64) -> f64 {
      let i = t.iter().position(|&ti| ti > x).unwrap_or(t.len() - 1).max(1);
      y[i - 1] + (x - t[i - 1]) / (t[i] - t[i - 1]) * (y[i] - y[i - 1])
    }
  }
}
