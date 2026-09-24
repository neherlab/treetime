#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionNegLog;
  use crate::distribution_core::function::DistributionFunction;
  use crate::{Distribution, convolve_across_edge};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, Side, SoftTailLaw};

  #[test]
  fn test_edge_convolution_variances_add_on_mass_window() -> Result<(), Report> {
    let result = convolve_across_edge(&gaussian(1.0)?, &gaussian(1.0)?, Side::Right, 1e-6, 400)?;

    let y0 = result.eval(0.0)?;
    assert_abs_diff_eq!(0.0, result.likely_time().unwrap(), epsilon = 0.05);
    assert_abs_diff_eq!(0.25, result.eval(1.0)? - y0, epsilon = 1e-6);
    assert_abs_diff_eq!(1.00, result.eval(2.0)? - y0, epsilon = 1e-5);
    Ok(())
  }

  #[test]
  fn test_edge_convolution_sets_soft_and_hard_tails() -> Result<(), Report> {
    for (soft, hard_side) in [(Side::Left, Side::Right), (Side::Right, Side::Left)] {
      let result = convolve_across_edge(&gaussian(1.0)?, &gaussian(1.0)?, soft, 1e-6, 200)?;
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

  fn gaussian(sigma: f64) -> Result<DistributionNegLog, Report> {
    let t = Array1::linspace(-6.0, 6.0, 241);
    let y = t.mapv(|t| t * t / (2.0 * sigma * sigma));
    let f = DistributionFunction::from_start_dx_values(t[0], t[1] - t[0], y)?;
    let left = SoftTailLaw::fit(f.grid_fn(), Side::Left, DEFAULT_TAIL_FIT_POINTS)?;
    let right = SoftTailLaw::fit(f.grid_fn(), Side::Right, DEFAULT_TAIL_FIT_POINTS)?;
    let f = f
      .with_left_extrap(BoundaryBehavior::Linear(left))?
      .with_right_extrap(BoundaryBehavior::Linear(right))?;
    Ok(Distribution::Function(f))
  }
}
