#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionNegLog;
  use crate::distribution_core::function::DistributionFunction;
  use crate::distribution_ops::convolve::distribution_convolution_fine;
  use crate::distribution_ops::mass_domain::peak_normalized_if_mass_sizable;
  use crate::policy::NegLog;
  use crate::{
    Distribution, MassProfile, TailMass, convolve_across_edge, convolve_forward_across_edge,
    distribution_multiplication,
  };
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_grid::MaxGridPoints;
  use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, HardApproachLaw, Side, SoftTailLaw};
  use treetime_utils::cast::f64_to_usize_saturating;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  #[test]
  fn test_edge_convolution_variances_add_on_mass_window() -> Result<(), Report> {
    let result = convolve_across_edge(
      &gaussian(1.0)?,
      &gaussian(1.0)?,
      Side::Right,
      1e-6,
      400,
      MaxGridPoints::default(),
    )?;

    let y0 = result.eval(0.0)?;
    assert_abs_diff_eq!(0.0, result.likely_time()?.unwrap(), epsilon = 0.05);
    assert_abs_diff_eq!(0.25, result.eval(1.0)? - y0, epsilon = 1e-6);
    assert_abs_diff_eq!(1.00, result.eval(2.0)? - y0, epsilon = 1e-5);
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
    )?;

    let Distribution::Function(f) = &result else {
      panic!("expected a gridded Function message");
    };
    assert_abs_diff_eq!(0.0, result.likely_time()?.unwrap(), epsilon = f.dx());
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

    let result = convolve_across_edge(&fine, &coarse, soft, 5e-4, 300, MaxGridPoints::default())?;

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
      distribution_convolution_fine(&narrow_fine, &wide_coarse, MaxGridPoints::default())?
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
  #[case::posterior_as_wide_as_the_message(2.0,  false)]
  #[case::posterior_narrower_than_the_cap( 0.05, true )]
  #[trace]
  fn test_edge_convolution_forward_message_resamples_to_resolve_a_narrow_posterior(
    #[case] posterior_sigma: f64,
    #[case] refined: bool,
  ) -> Result<(), Report> {
    let forward = convolve_forward_across_edge(
      &gaussian_at(1.0, 0.001)?,
      &gaussian_at(1.0, 0.02)?,
      5e-4,
      300,
      MaxGridPoints::default(),
    )?;
    let capped = forward.capped()?;
    let posterior = gaussian_at(posterior_sigma, 0.001)?;

    let resolving = forward.resampled_to_resolve(&posterior, 1e-9)?;

    let Distribution::Function(capped) = &capped else {
      panic!("expected a gridded Function message");
    };
    let window = capped.x_max() - capped.x_min();
    let intervals = (window / 0.02).ceil();
    assert!(window / 299.0 > 0.02, "the window {window} needs more than 300 points at the cap");
    assert_eq!(f64_to_usize_saturating(intervals) + 1, capped.len());
    pretty_assert_abs_diff_eq!(window / intervals, capped.dx(), epsilon = 1e-12);
    let Distribution::Function(p) = &posterior else {
      panic!("expected a gridded Function posterior");
    };
    let profile = MassProfile::new(p)?;
    let fractions = profile.fractions(TailMass::Required).expect("finite mass");
    let resolving_dx = (fractions.survival_inverse(1e-9) - fractions.quantile(1e-9)) / 299.0;
    match resolving {
      None => {
        assert!(!refined, "a posterior {posterior_sigma} wide must take a finer message");
        assert!(resolving_dx >= capped.dx());
      },
      Some(Distribution::Function(f)) => {
        assert!(refined, "a posterior {posterior_sigma} wide keeps the capped message");
        let width = f.x_max() - f.x_min();
        pretty_assert_abs_diff_eq!(capped.x_min(), f.x_min(), epsilon = 1e-12);
        pretty_assert_abs_diff_eq!(capped.x_max(), f.x_max(), epsilon = 1e-12);
        assert_eq!(f64_to_usize_saturating((width / resolving_dx).ceil()) + 1, f.len());
      },
      Some(other) => panic!("expected a gridded Function message, got {other:?}"),
    }
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
