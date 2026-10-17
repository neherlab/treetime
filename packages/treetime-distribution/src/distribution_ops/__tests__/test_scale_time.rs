#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionNegLog;
  use crate::DistributionFunction;
  use crate::policy::NegLog;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use treetime_grid::{BoundaryBehavior, HardApproachLaw, SoftTailLaw};
  use treetime_utils::assert_error;

  #[test]
  fn test_scale_time_moves_values_and_keeps_level() -> Result<(), Report> {
    let f = DistributionNegLog::Function(helpers::branch_like()?.with_amplitude(1.5));
    let scaled = f.scale_time(0.25)?;
    for t in [0.05, 0.3, 1.0, 2.5, 3.9, 6.0] {
      assert_abs_diff_eq!(f.eval(t)?, scaled.eval(0.25 * t)?, epsilon = 1e-12);
    }
    Ok(())
  }

  #[test]
  fn test_scale_time_scales_the_mass() -> Result<(), Report> {
    let f = DistributionNegLog::Function(helpers::branch_like()?);
    assert_abs_diff_eq!(
      f.log_mass()? + 3.0_f64.ln(),
      f.scale_time(3.0)?.log_mass()?,
      epsilon = 1e-12
    );
    Ok(())
  }

  #[test]
  fn test_scale_time_transforms_the_tail_laws() -> Result<(), Report> {
    let DistributionNegLog::Function(scaled) = DistributionNegLog::Function(helpers::branch_like()?).scale_time(2.0)?
    else {
      panic!("expected a function");
    };
    assert_eq!(
      BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 2.0 }),
      scaled.left_extrap()
    );
    assert_eq!(
      BoundaryBehavior::Linear(SoftTailLaw { slope: 0.75 }),
      scaled.right_extrap()
    );
    Ok(())
  }

  #[test]
  fn test_scale_time_point_range_and_invalid_factor() -> Result<(), Report> {
    assert_eq!(
      DistributionNegLog::point(4.0, 0.5),
      DistributionNegLog::point(2.0, 0.5).scale_time(2.0)?
    );
    assert_eq!(
      DistributionNegLog::range((2.0, 6.0), 0.5),
      DistributionNegLog::range((1.0, 3.0), 0.5).scale_time(2.0)?
    );
    assert_error!(
      DistributionNegLog::point(2.0, 0.5).scale_time(0.0),
      "Time scale factor must be positive, got 0"
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    /// Branch-like likelihood t^2 exp(-1.5 t) on [0.1, 4]: power-law approach to t = 0, exponential right tail
    pub(super) fn branch_like() -> Result<DistributionFunction<f64, NegLog>, Report> {
      let t = Array1::linspace(0.1, 4.0, 391);
      let y = t.mapv(|t: f64| 1.5 * t - 2.0 * t.ln());
      DistributionFunction::from_start_dx_values(t[0], t[1] - t[0], y)?
        .with_left_extrap(BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 2.0 }))?
        .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 1.5 }))
    }
  }
}
