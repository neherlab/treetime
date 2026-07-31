#[cfg(test)]
mod tests {
  use crate::distribution_core::function::{DistributionFunction, TailSide};
  use crate::policy::Plain;
  use eyre::Report;
  use ndarray::Array1;
  use treetime_utils::pretty_assert_ulps_eq;

  type DistFn = DistributionFunction<f64, Plain>;

  const A: f64 = 5.0;
  const DX: f64 = 0.1;

  /// A pure exponential `exp(-A x)` sampled on `[0, 1]` at spacing `DX`. Its logarithm is exactly
  /// linear, so the boundary-slope extrapolation reproduces the analytic curve with no model error,
  /// making it the oracle for the decay math.
  fn exponential_grid() -> Result<DistFn, Report> {
    let y = Array1::from_shape_fn(11, |i| (-A * DX * i as f64).exp());
    DistributionFunction::from_start_dx_values(0.0, DX, y)
  }

  #[test]
  fn test_tail_extension_right_matches_analytic_exponential() -> Result<(), Report> {
    let f = exponential_grid()?;
    let extended = f.extend_log_linear_tail(TailSide::Right, 1e-10, 10_000)?;

    // The original grid is preserved as a prefix and the appended points continue exp(-A x).
    let x = extended.t();
    let expected = x.mapv(|xi| (-A * xi).exp());
    pretty_assert_ulps_eq!(extended.y().clone(), expected, max_ulps = 64);

    // The grid grew only on the right: x_min unchanged, x_max pushed out.
    pretty_assert_ulps_eq!(extended.x_min(), 0.0);
    assert!(extended.x_max() > f.x_max());
    Ok(())
  }

  #[test]
  fn test_tail_extension_right_terminates_at_relative_floor() -> Result<(), Report> {
    let rel_floor = 1e-10;
    let f = exponential_grid()?;
    let extended = f.extend_log_linear_tail(TailSide::Right, rel_floor, 10_000)?;

    let peak = f.y().iter().copied().fold(f64::MIN, f64::max);
    let floor = rel_floor * peak;

    // Every retained point is at or above the floor, and the next point would fall below it:
    // extension stopped because the density became negligible, not because of a budget cap.
    let boundary = *extended.y().last().unwrap();
    assert!(boundary >= floor, "boundary {boundary:e} below floor {floor:e}");
    let next = boundary * (-A * DX).exp();
    assert!(next < floor, "next point {next:e} not below floor {floor:e}");
    Ok(())
  }

  #[test]
  fn test_tail_extension_flat_boundary_is_unchanged() -> Result<(), Report> {
    // A constant density does not decay outward, so neither side is extended.
    let f: DistFn = DistributionFunction::from_start_dx_values(0.0, DX, Array1::from_elem(11, 2.0))?;
    let extended = f.extend_log_linear_tail(TailSide::Right, 1e-10, 10_000)?;
    pretty_assert_ulps_eq!(extended.y().clone(), f.y().clone());
    assert_eq!(extended.len(), f.len());
    Ok(())
  }

  #[test]
  fn test_tail_extension_rising_boundary_is_unchanged() -> Result<(), Report> {
    // Extending the side where the density rises outward (here, the left of a decreasing curve)
    // is a no-op: the boundary value exceeds the interior anchor, so there is nothing to decay.
    let f = exponential_grid()?;
    let extended = f.extend_log_linear_tail(TailSide::Left, 1e-10, 10_000)?;
    pretty_assert_ulps_eq!(extended.y().clone(), f.y().clone());
    assert_eq!(extended.len(), f.len());
    Ok(())
  }

  #[test]
  fn test_tail_extension_right_leaves_left_side_untouched() -> Result<(), Report> {
    let f = exponential_grid()?;
    let extended = f.extend_log_linear_tail(TailSide::Right, 1e-10, 10_000)?;
    // The original samples remain the leading prefix, bit-for-bit.
    pretty_assert_ulps_eq!(extended.y().slice(ndarray::s![..f.len()]).to_owned(), f.y().clone());
    pretty_assert_ulps_eq!(extended.x_min(), f.x_min());
    Ok(())
  }

  #[test]
  fn test_tail_extension_left_matches_analytic_and_leaves_right_untouched() -> Result<(), Report> {
    // A left-decaying curve exp(A x) on [0, 1]: decays as x decreases, so Left extends and Right
    // (rising outward) is a no-op.
    let y = Array1::from_shape_fn(11, |i| (A * DX * i as f64).exp());
    let f: DistFn = DistributionFunction::from_start_dx_values(0.0, DX, y)?;

    let extended = f.extend_log_linear_tail(TailSide::Left, 1e-10, 10_000)?;
    let x = extended.t();
    let expected = x.mapv(|xi| (A * xi).exp());
    pretty_assert_ulps_eq!(extended.y().clone(), expected, max_ulps = 64);

    assert!(extended.x_min() < f.x_min());
    pretty_assert_ulps_eq!(extended.x_max(), f.x_max());
    Ok(())
  }

  #[test]
  fn test_tail_extension_respects_point_budget() -> Result<(), Report> {
    let f = exponential_grid()?;
    let budget = 5;
    let extended = f.extend_log_linear_tail(TailSide::Right, 1e-10, budget)?;
    // The budget caps how many points are appended, short of the floor.
    assert_eq!(extended.len(), f.len() + budget);
    Ok(())
  }

  #[test]
  fn test_tail_extension_zero_budget_is_unchanged() -> Result<(), Report> {
    let f = exponential_grid()?;
    let extended = f.extend_log_linear_tail(TailSide::Right, 1e-10, 0)?;
    assert_eq!(extended.len(), f.len());
    Ok(())
  }

  #[test]
  fn test_tail_extension_floor_above_boundary_is_unchanged() -> Result<(), Report> {
    // A floor at the peak leaves the boundary already below it, so nothing is appended.
    let f = exponential_grid()?;
    let extended = f.extend_log_linear_tail(TailSide::Right, 1.0, 10_000)?;
    assert_eq!(extended.len(), f.len());
    Ok(())
  }
}
