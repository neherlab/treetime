#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionNegLog;
  use crate::distribution_core::function::DistributionFunction;
  use crate::distribution_ops::convolve::distribution_convolution_fine;
  use crate::distribution_ops::mass_domain::{peak_normalized_if_mass_sizable, trim_to_mass};
  use crate::policy::NegLog;
  use crate::{
    Distribution, MassProfile, MessageSpacing, TailMass, convolve_across_edge, convolve_forward_across_edge,
    distribution_multiplication,
  };
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::{Array1, Axis};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use statrs::function::erf::{erf, erfc};
  use std::f64::consts::{LN_10, PI, SQRT_2};
  use treetime_grid::MaxGridPoints;
  use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, HardApproachLaw, Side, SoftTailLaw};
  use treetime_ops::FftPlans;
  use treetime_utils::array::ndarray::min_or;
  use treetime_utils::cast::f64_to_usize_saturating;
  use treetime_utils::{make_error, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  #[test]
  fn test_edge_convolution_variances_add_on_mass_window() -> Result<(), Report> {
    let result = convolve_across_edge(
      &gaussian(1.0)?,
      &gaussian(1.0)?,
      Side::Right,
      1e-6,
      400,
      MaxGridPoints::default(),
      &FftPlans::new(),
    )?;

    let y0 = result.eval(0.0)?;
    pretty_assert_abs_diff_eq!(0.0, result.likely_time()?.unwrap(), epsilon = 0.05);
    pretty_assert_abs_diff_eq!(0.25, result.eval(1.0)? - y0, epsilon = 1e-6);
    pretty_assert_abs_diff_eq!(1.00, result.eval(2.0)? - y0, epsilon = 1e-5);
    Ok(())
  }

  #[test]
  fn test_edge_convolution_sets_soft_and_hard_tails() -> Result<(), Report> {
    for (soft, hard_side) in [(Side::Left, Side::Right), (Side::Right, Side::Left)] {
      let result = convolve_across_edge(
        &gaussian(1.0)?,
        &gaussian(1.0)?,
        soft,
        1e-6,
        200,
        MaxGridPoints::default(),
        &FftPlans::new(),
      )?;
      let Distribution::Function(f) = &result else {
        panic!("expected a gridded Function message");
      };
      assert!(f.len() >= 200, "output holds at least grid_points points");
      let (soft_extrap, hard_extrap) = match soft {
        Side::Left => (f.left_extrap(), f.right_extrap()),
        Side::Right => (f.right_extrap(), f.left_extrap()),
      };
      assert!(
        matches!(soft_extrap, BoundaryBehavior::Linear(_)),
        "soft {soft:?} side is a fitted log-linear tail"
      );
      assert_eq!(
        BoundaryBehavior::Hard,
        hard_extrap,
        "hard {hard_side:?} side is the Minkowski bound"
      );
    }
    Ok(())
  }

  #[test]
  fn test_edge_convolution_operand_without_mass_uses_result_window() -> Result<(), Report> {
    let untailed = Distribution::Function(gaussian_grid(1.0)?);
    assert!(peak_normalized_if_mass_sizable(&gaussian_grid(1.0)?).is_none());

    let result = convolve_across_edge(
      &untailed,
      &gaussian(1.0)?,
      Side::Right,
      1e-6,
      400,
      MaxGridPoints::default(),
      &FftPlans::new(),
    )?;

    let Distribution::Function(f) = &result else {
      panic!("expected a gridded Function message");
    };
    pretty_assert_abs_diff_eq!(0.0, result.likely_time()?.unwrap(), epsilon = f.dx());
    Ok(())
  }

  #[test]
  fn test_edge_convolution_point_keeps_steep_soft_law_without_window() -> Result<(), Report> {
    let steep = DistributionFunction::<f64, NegLog>::from_range_values((0.0, 1.0), Array1::zeros(2))?
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 1e6 }))?;

    let message = convolve_across_edge(
      &DistributionNegLog::point(0.0, 0.0),
      &Distribution::Function(steep.clone()),
      Side::Right,
      5e-4,
      300,
      MaxGridPoints::default(),
      &FftPlans::new(),
    )?;

    assert_eq!(Distribution::Function(steep), message);
    Ok(())
  }

  #[test]
  fn test_edge_convolution_point_keeps_the_level() -> Result<(), Report> {
    let branch = branch(BoundaryBehavior::Hard, 0.0)?;

    let message = convolve_across_edge(
      &DistributionNegLog::point(2000.0, 1.5),
      &branch.negate()?,
      Side::Left,
      5e-4,
      300,
      MaxGridPoints::default(),
      &FftPlans::new(),
    )?;

    pretty_assert_abs_diff_eq!(branch.log_mass()? - 1.5, message.log_mass()?, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_edge_convolution_date_range_message_ends_at_range_end() -> Result<(), Report> {
    let branch = branch(BoundaryBehavior::Hard, 0.0)?;
    let leaf = DistributionNegLog::range((2000.0, 2001.0), 0.0);

    let message = convolve_across_edge(
      &leaf,
      &branch.negate()?,
      Side::Left,
      1e-6,
      300,
      MaxGridPoints::default(),
      &FftPlans::new(),
    )?;

    let Distribution::Function(f) = &message else {
      panic!("expected a gridded Function message, got {message:?}");
    };
    pretty_assert_ulps_eq!(2001.0, f.x_max(), max_ulps = 0);
    assert_eq!(BoundaryBehavior::Hard, f.right_extrap());
    let inside_last_cell = (-message.eval(2001.0 - f.dx() / 2.0)?).exp();
    assert!(
      inside_last_cell > 0.0,
      "density {inside_last_cell} inside the last cell"
    );
    Ok(())
  }

  #[test]
  fn test_edge_convolution_exact_dated_leaf_messages_multiply() -> Result<(), Report> {
    let negated = branch(
      BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 2.0 }),
      0.1,
    )?
    .negate()?;
    let messages = [2000.0, 2000.5].map(|date| {
      convolve_across_edge(
        &DistributionNegLog::point(date, 0.0),
        &negated,
        Side::Left,
        1e-6,
        300,
        MaxGridPoints::default(),
        &FftPlans::new(),
      )
    });
    let [Ok(first), Ok(second)] = messages else {
      panic!("expected two messages, got {messages:?}");
    };

    for message in [&first, &second] {
      let Distribution::Function(f) = message else {
        panic!("expected a gridded Function message, got {message:?}");
      };
      assert_eq!(BoundaryBehavior::Hard, f.right_extrap());
    }
    let product = distribution_multiplication(&first, &second, MaxGridPoints::default())?;
    assert!(
      matches!(product, Distribution::Function(_)),
      "expected a gridded product, got {product:?}"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::forward_window_wider_than_the_cap(  Side::Right, 1.0, 0.01 )]
  #[case::forward_window_narrower_than_the_cap(Side::Right, 0.1, 0.01 )]
  #[case::backward_keeps_the_finer_spacing(    Side::Left,  1.0, 0.001)]
  #[trace]
  fn test_edge_convolution_output_spacing_is_at_most_the_spacing_cap(
    #[case] soft: Side,
    #[case] sigma: f64,
    #[case] max_dx: f64,
  ) -> Result<(), Report> {
    let fine = gaussian_at(sigma, 0.001)?;
    let coarse = gaussian_at(sigma, 0.01)?;

    let result = convolve_across_edge(&fine, &coarse, soft, 5e-4, 300, MaxGridPoints::default(), &FftPlans::new())?;

    let Distribution::Function(f) = &result else {
      panic!("expected a gridded Function message");
    };
    let width = f.x_max() - f.x_min();
    let expected = if width / 299.0 <= max_dx {
      300
    } else {
      f64_to_usize_saturating((width / max_dx).ceil()) + 1
    };
    assert_eq!(expected, f.len());
    Ok(())
  }

  #[test]
  fn test_edge_convolution_forward_tail_is_the_gaussian_slope_at_the_mean_fit_position() -> Result<(), Report> {
    let narrow_fine = gaussian_at(0.1, 0.001)?;
    let wide_coarse = gaussian_at(1.0, 0.02)?;
    let Distribution::Function(fine) =
      distribution_convolution_fine(&narrow_fine, &wide_coarse, MaxGridPoints::default(), &FftPlans::new())?
    else {
      panic!("expected a gridded convolution");
    };

    let result = convolve_across_edge(
      &narrow_fine,
      &wide_coarse,
      Side::Right,
      5e-4,
      300,
      MaxGridPoints::default(),
      &FftPlans::new(),
    )?;

    let Distribution::Function(f) = &result else {
      panic!("expected a gridded Function message");
    };
    let width = f.x_max() - f.x_min();
    assert!(
      width / 299.0 > 0.02,
      "the window of {width} needs more than 300 points at the cap"
    );
    assert_eq!(f64_to_usize_saturating((width / 0.02).ceil()) + 1, f.len());
    let law = f.right_extrap().soft_law().expect("a soft right tail");
    let t = fine.t();
    let edge = t
      .iter()
      .rposition(|&x| x <= f.x_max())
      .expect("a convolution point below the window edge");
    let variance = 0.1 * 0.1 + 1.0 * 1.0;
    pretty_assert_abs_diff_eq!((t[edge] - 2.0 * fine.dx()) / variance, law.slope, epsilon = 1e-8);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::posterior_as_wide_as_the_message(2.0,  MessageSpacing::Capped   )]
  #[case::posterior_narrower_than_the_cap( 0.05, MessageSpacing::Resolving)]
  #[trace]
  fn test_edge_convolution_forward_message_resamples_to_resolve_a_narrow_posterior(
    #[case] posterior_sigma: f64,
    #[case] expected_spacing: MessageSpacing,
  ) -> Result<(), Report> {
    let forward = convolve_forward_across_edge(
      &gaussian_at(1.0, 0.001)?,
      &gaussian_at(1.0, 0.02)?,
      5e-4,
      300,
      MaxGridPoints::default(), &FftPlans::new(),
    )?;
    let posterior = gaussian_at(posterior_sigma, 0.001)?;

    let resolved = forward
      .resolve_posterior(1e-9, |_| Ok(posterior.clone()))?
      .into_untrimmed();

    let Distribution::Function(p) = &posterior else {
      panic!("expected a gridded Function posterior");
    };
    let profile = MassProfile::new(p)?;
    let fractions = profile.fractions(TailMass::Required).expect("finite mass");
    let resolving_dx = (fractions.survival_inverse(1e-9) - fractions.quantile(1e-9)) / 299.0;
    let max_dx = match expected_spacing {
      MessageSpacing::Capped => 0.02,
      MessageSpacing::Exact | MessageSpacing::Resolving | MessageSpacing::Fine => resolving_dx,
    };
    let Distribution::Function(f) = &resolved.message else {
      panic!("expected a gridded Function message");
    };
    let window = f.x_max() - f.x_min();
    let intervals = (window / max_dx).ceil();
    assert_eq!(expected_spacing, resolved.spacing);
    assert!(window / 299.0 > max_dx, "the window {window} needs more than 300 points");
    assert_eq!(f64_to_usize_saturating(intervals) + 1, f.len());
    pretty_assert_abs_diff_eq!(window / intervals, f.dx(), epsilon = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::posterior_as_wide_as_the_message(2.0,  MessageSpacing::Capped   )]
  #[case::posterior_narrower_than_the_cap( 0.05, MessageSpacing::Resolving)]
  #[trace]
  fn test_edge_convolution_forward_posterior_trims_like_trim_to_mass(
    #[case] posterior_sigma: f64,
    #[case] expected_spacing: MessageSpacing,
  ) -> Result<(), Report> {
    let forward = convolve_forward_across_edge(
      &gaussian_at(1.0, 0.001)?,
      &gaussian_at(1.0, 0.02)?,
      5e-4,
      300,
      MaxGridPoints::default(), &FftPlans::new(),
    )?;
    let posterior = gaussian_at(posterior_sigma, 0.001)?;

    let trimmed = forward
      .resolve_posterior(1e-6, |_| Ok(posterior.clone()))?
      .trim_posterior_to_mass(1e-6)?;

    let expected = trim_to_mass(posterior.clone(), 1e-6)?;
    assert_eq!(expected_spacing, trimmed.spacing);
    assert!(expected != posterior, "the posterior must lose mass beyond its 1e-6 quantiles");
    assert_eq!(expected, trimmed.posterior);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::backward(Side::Left)]
  #[case::forward( Side::Right)]
  #[trace]
  fn test_edge_convolution_bimodal_operand_with_deep_valley_keeps_both_modes(#[case] soft: Side) -> Result<(), Report> {
    let result = convolve_across_edge(
      &two_gaussians(20.0, 80.0, 1.0, 1001)?,
      &gaussian(1.0)?,
      soft,
      1e-6,
      400,
      MaxGridPoints::default(), &FftPlans::new(),
    )?;

    let y20 = result.eval(20.0)?;
    pretty_assert_abs_diff_eq!(0.0, result.eval(80.0)? - y20, epsilon = 1e-6);
    let valley = result.eval(50.0)? - y20;
    assert!(
      valley.is_finite() && valley > 20.0,
      "the valley at t = 50 must stay finite and more than 20 nats below the modes, got {valley}"
    );
    Ok(())
  }

  #[test]
  fn test_convolution_of_bimodal_operand_with_deep_valley_matches_its_oracles() -> Result<(), Report> {
    let (a, b) = (two_gaussians(20.0, 80.0, 1.0, 1441)?, gaussian(1.0)?);
    let conv = distribution_convolution_fine(&a, &b, MaxGridPoints::default(), &FftPlans::new())?;
    let Distribution::Function(conv) = conv else {
      panic!("the convolution of two grid functions must be a grid function, got {conv}");
    };
    let t = conv.t();
    let y = conv.y();
    assert!(y.iter().all(|y| !y.is_nan()), "the convolution holds NaN: {y}");
    let exact = t.mapv(deep_valley_exact_neg_log);
    let quadrature = log_space_trapezoid_convolution(&a, &b, &t)?;
    let exact_peak = min_or(&exact, f64::INFINITY);

    let mode_samples = t
      .iter()
      .positions(|&t| (t - 20.0).abs() < 1e-9 || (t - 80.0).abs() < 1e-9)
      .collect_vec();
    assert_eq!(2, mode_samples.len());
    for i in mode_samples {
      pretty_assert_abs_diff_eq!(exact[i], y[i], epsilon = 1e-10);
    }

    let neg_log_total = -(conv.dx() * y.mapv(|y| (-y).exp()).sum()).ln();
    pretty_assert_abs_diff_eq!(deep_valley_exact_neg_log_total(), neg_log_total, epsilon = 1e-10);

    let near_modes = (0..t.len())
      .filter(|&i| (t[i] - 20.0).abs() <= 6.0 || (t[i] - 80.0).abs() <= 6.0)
      .collect_vec();
    pretty_assert_abs_diff_eq!(
      exact.select(Axis(0), &near_modes),
      y.select(Axis(0), &near_modes),
      epsilon = 1e-6
    );

    let within_fft_range = (0..t.len())
      .filter(|&i| quadrature[i] - exact_peak <= 20.0)
      .collect_vec();
    pretty_assert_abs_diff_eq!(
      quadrature.select(Axis(0), &within_fft_range),
      y.select(Axis(0), &within_fft_range),
      epsilon = 1e-6
    );

    let below_fft_range = (0..t.len())
      .filter(|&i| t[i] > 20.0 && t[i] < 80.0 && quadrature[i] - exact_peak > CONV_TRUST_DEPTH + 1.0)
      .collect_vec();
    assert!(!below_fft_range.is_empty());
    pretty_assert_abs_diff_eq!(
      quadrature.select(Axis(0), &below_fft_range),
      y.select(Axis(0), &below_fft_range),
      epsilon = 1e-9
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn gaussian_at(sigma: f64, dx: f64) -> Result<DistributionNegLog, Report> {
      let n_points = f64_to_usize_saturating((12.0 * sigma / dx).round()) + 1;
      let t = Array1::linspace(-6.0 * sigma, 6.0 * sigma, n_points);
      let y = t.mapv(|t| t * t / (2.0 * sigma * sigma));
      let f = DistributionFunction::from_start_dx_values(t[0], t[1] - t[0], y)?;
      let left = SoftTailLaw::fit(f.grid_fn(), Side::Left, DEFAULT_TAIL_FIT_POINTS)?;
      let right = SoftTailLaw::fit(f.grid_fn(), Side::Right, DEFAULT_TAIL_FIT_POINTS)?;
      let f = f
        .with_left_extrap(BoundaryBehavior::Linear(left))?
        .with_right_extrap(BoundaryBehavior::Linear(right))?;
      Ok(Distribution::Function(f))
    }

    pub(super) fn two_gaussians(
      mu_a: f64,
      mu_b: f64,
      sigma: f64,
      n_points: usize,
    ) -> Result<DistributionNegLog, Report> {
      let t = Array1::linspace(mu_a - 6.0 * sigma, mu_b + 6.0 * sigma, n_points);
      let y = t.mapv(|t| ((t - mu_a).powi(2)).min((t - mu_b).powi(2)) / (2.0 * sigma * sigma));
      let f = DistributionFunction::from_start_dx_values(t[0], t[1] - t[0], y)?;
      let left = SoftTailLaw::fit(f.grid_fn(), Side::Left, DEFAULT_TAIL_FIT_POINTS)?;
      let right = SoftTailLaw::fit(f.grid_fn(), Side::Right, DEFAULT_TAIL_FIT_POINTS)?;
      let f = f
        .with_left_extrap(BoundaryBehavior::Linear(left))?
        .with_right_extrap(BoundaryBehavior::Linear(right))?;
      Ok(Distribution::Function(f))
    }

    pub(super) const CONV_TRUST_DEPTH: f64 = 13.0 * LN_10;

    pub(super) fn log_space_trapezoid_convolution(
      a: &DistributionNegLog,
      b: &DistributionNegLog,
      t: &Array1<f64>,
    ) -> Result<Array1<f64>, Report> {
      let (Distribution::Function(a), Distribution::Function(b)) = (a, b) else {
        return make_error!("the quadrature oracle needs two grid functions");
      };
      let dx = a.dx();
      let (ya, yb) = (a.y(), b.y());
      let (last_a, last_b) = (ya.len() - 1, yb.len() - 1);
      let t0 = a.x_min() + b.x_min();
      t.iter()
        .map(|&t| {
          let k = f64_to_usize_saturating(((t - t0) / dx).round());
          let (lo, hi) = (k.saturating_sub(last_b), k.min(last_a));
          let terms = (lo..=hi)
            .map(|i| {
              let weight: f64 = if lo == hi {
                0.25
              } else if i == lo || i == hi {
                0.5
              } else {
                1.0
              };
              weight.ln() - ya[i] - yb[k - i]
            })
            .collect_vec();
          let max = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
          Ok(-(max + terms.iter().map(|v| (v - max).exp()).sum::<f64>().ln() + dx.ln()))
        })
        .collect()
    }

    pub(super) fn deep_valley_exact_neg_log(t: f64) -> f64 {
      let log_masses = [(20.0, 14.0, 50.0), (80.0, 50.0, 86.0)].map(|(mu, s_lo, s_hi)| {
        let d = t - mu;
        let (u_lo, u_hi) = ((t - s_hi).max(-6.0), (t - s_lo).min(6.0));
        -d * d / 4.0 + ln_gauss_integral(u_lo - d / 2.0, u_hi - d / 2.0)
      });
      let [a, b] = log_masses;
      let (hi, lo) = (a.max(b), a.min(b));
      if hi == f64::NEG_INFINITY {
        return f64::INFINITY;
      }
      -(hi + (lo - hi).exp().ln_1p())
    }

    pub(super) fn deep_valley_exact_neg_log_total() -> f64 {
      let half_sqrt_2pi = (PI / 2.0).sqrt();
      let mode_mass = half_sqrt_2pi * (erf(30.0 / SQRT_2) + erf(6.0 / SQRT_2));
      let kernel_mass = half_sqrt_2pi * 2.0 * erf(6.0 / SQRT_2);
      -(2.0 * mode_mass * kernel_mass).ln()
    }

    fn ln_gauss_integral(lo: f64, hi: f64) -> f64 {
      let ln_half_sqrt_pi = (PI.sqrt() / 2.0).ln();
      if hi <= lo {
        f64::NEG_INFINITY
      } else if lo >= 0.0 {
        ln_half_sqrt_pi + erfc(lo).ln() + (-erfc(hi) / erfc(lo)).ln_1p()
      } else if hi <= 0.0 {
        ln_gauss_integral(-hi, -lo)
      } else {
        ln_half_sqrt_pi + (erf(hi) - erf(lo)).ln()
      }
    }

    pub(super) fn gaussian(sigma: f64) -> Result<DistributionNegLog, Report> {
      let f = gaussian_grid(sigma)?;
      let left = SoftTailLaw::fit(f.grid_fn(), Side::Left, DEFAULT_TAIL_FIT_POINTS)?;
      let right = SoftTailLaw::fit(f.grid_fn(), Side::Right, DEFAULT_TAIL_FIT_POINTS)?;
      let f = f
        .with_left_extrap(BoundaryBehavior::Linear(left))?
        .with_right_extrap(BoundaryBehavior::Linear(right))?;
      Ok(Distribution::Function(f))
    }

    pub(super) fn branch(left: BoundaryBehavior, x_min: f64) -> Result<DistributionNegLog, Report> {
      let t = Array1::linspace(x_min, x_min + 2.0, 41);
      let f = DistributionFunction::<f64, NegLog>::from_range_values(
        (x_min, x_min + 2.0),
        t.mapv(|t| 4.0 * (t - 0.5).powi(2)),
      )?
      .with_left_extrap(left)?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 6.0 }))?;
      Ok(Distribution::Function(f))
    }

    pub(super) fn gaussian_grid(sigma: f64) -> Result<DistributionFunction<f64, NegLog>, Report> {
      let t = Array1::linspace(-6.0, 6.0, 241);
      let y = t.mapv(|t| t * t / (2.0 * sigma * sigma));
      DistributionFunction::from_start_dx_values(t[0], t[1] - t[0], y)
    }
  }

  use helpers::*;
}
