#[cfg(test)]
mod tests {
  use crate::*;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::at_x_min((0.0, 1.0), array![0.0, 10.0], 0.0, 0.0)]
  #[case::at_x_max((0.0, 1.0), array![0.0, 10.0], 1.0, 10.0)]
  #[case::at_midpoint((0.0, 2.0), array![0.0, 5.0, 10.0], 1.0, 5.0)]
  #[case::negative_range_at_min((-5.0, 5.0), array![100.0, 200.0], -5.0, 100.0)]
  #[case::negative_range_at_max((-5.0, 5.0), array![100.0, 200.0], 5.0, 200.0)]
  #[trace]
  fn test_gridfn_interp_exact_grid_points(
    #[case] x_range: (f64, f64),
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, max_ulps = 4);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::midpoint_linear((0.0, 1.0), array![0.0, 10.0], 0.5, 5.0)]
  #[case::at_30_percent((0.0, 1.0), array![0.0, 10.0], 0.3, 3.0)]
  #[case::at_70_percent((0.0, 1.0), array![0.0, 10.0], 0.7, 7.0)]
  #[case::midpoint_offset((0.0, 2.0), array![10.0, 20.0], 1.0, 15.0)]
  #[case::negative_to_positive((-1.0, 1.0), array![0.0, 100.0], 0.0, 50.0)]
  #[case::three_points_first_half((0.0, 2.0), array![0.0, 10.0, 20.0], 0.5, 5.0)]
  #[case::three_points_second_half((0.0, 2.0), array![0.0, 10.0, 20.0], 1.5, 15.0)]
  #[trace]
  fn test_gridfn_interp_interior_points(
    #[case] x_range: (f64, f64),
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, max_ulps = 4);
    Ok(())
  }

  // `resample_range_dx_clamped` holds the boundary value where grid-construction rounding pushes the
  // final target point beyond the source support. `from_range_dx(0, 1, 0.4)` rounds to four points
  // ending at 1.2, past x_max = 1.0. A plain resample rejects that overshoot under the default `Error`
  // tail; the clamped resample clamps the query to x_max and reads the boundary ordinate there.
  #[test]
  fn test_gridfn_resample_range_dx_clamped_holds_boundary_on_overshoot() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, 10.0])?;
    assert_error!(
      grid_fn.resample_range_dx((0.0, 1.0), 0.4),
      "GridFn evaluated at 1.2000000000000002, above the support boundary 1.0, but no extrapolation policy is set for that side"
    );
    let clamped = grid_fn.resample_range_dx_clamped((0.0, 1.0), 0.4)?;
    // Target points 0.0, 0.4, 0.8, 1.2; the last clamps to x_max = 1.0, reading y = 10.0.
    assert_ulps_eq!(clamped.y(), &array![0.0, 4.0, 8.0, 10.0], max_ulps = 4);
    Ok(())
  }

  // Within the source support the clamp never fires, so a clamped resample matches a plain one.
  #[test]
  fn test_gridfn_resample_range_dx_clamped_matches_plain_within_support() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, 10.0])?;
    let plain = grid_fn.resample_range_dx((0.0, 1.0), 0.5)?;
    let clamped = grid_fn.resample_range_dx_clamped((0.0, 1.0), 0.5)?;
    assert_ulps_eq!(clamped.y(), plain.y(), max_ulps = 4);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  // Hard extrapolation: return 0.0 on either side (explicit opt-in)
  #[case::left((0.0, 1.0), array![2.0, 10.0], -1.0)]
  #[case::right((0.0, 1.0), array![2.0, 10.0], 2.0)]
  #[trace]
  fn test_gridfn_interp_extrapolation_zero(
    #[case] x_range: (f64, f64),
    #[case] y: Array1<f64>,
    #[case] query: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    let actual = grid_fn.interp_with_extrap(query, BoundaryBehavior::Hard, BoundaryBehavior::Hard)?;
    assert_ulps_eq!(0.0, actual, max_ulps = 4);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  // Default policy is Error: out-of-support evaluation fails on both sides.
  #[case::left(-1.0, "GridFn evaluated at -1.0, below the support boundary 0.0, but no extrapolation policy is set for that side")]
  #[case::right(2.0, "GridFn evaluated at 2.0, above the support boundary 1.0, but no extrapolation policy is set for that side")]
  #[trace]
  fn test_gridfn_interp_out_of_support_errors_by_default(#[case] query: f64, #[case] expected: &str) -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![2.0, 10.0])?;
    assert_error!(grid_fn.interp(query), expected);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_many() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![0.0, 10.0, 20.0])?;
    let queries = array![0.5, 1.0, 1.5];
    let expected = array![5.0, 10.0, 15.0];
    let actual = grid_fn.interp_many(&queries)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_from_grid() -> Result<(), Report> {
    let grid_fn = GridFn::from_grid((0.0, 1.0), 0.25, |x| x * x)?;
    assert_eq!(grid_fn.x().len(), 5);
    assert_ulps_eq!(grid_fn.x()[0], 0.0, max_ulps = 4);
    assert_ulps_eq!(grid_fn.x()[4], 1.0, max_ulps = 4);
    assert_ulps_eq!(grid_fn.y()[0], 0.0, max_ulps = 4);
    assert_ulps_eq!(grid_fn.y()[2], 0.25, max_ulps = 4);
    assert_ulps_eq!(grid_fn.y()[4], 1.0, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_gridfn_accessors() -> Result<(), Report> {
    let x_range = (0.0, 2.0);
    let y = array![10.0, 20.0, 30.0];
    let grid_fn = GridFn::from_range_values(x_range, y.clone())?;
    let x = grid_fn.x();
    assert_ulps_eq!(x[0], 0.0);
    assert_ulps_eq!(x[1], 1.0);
    assert_ulps_eq!(x[2], 2.0);
    assert_eq!(&y, grid_fn.y());
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    let expected = GridFn::from_range_values((0.0, 2.0), array![2.0, 4.0, 6.0])?;
    let actual = grid_fn.mapv(|y| y * 2.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_negate_arg_inplace() -> Result<(), Report> {
    let mut grid_fn = GridFn::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    grid_fn.negate_arg_inplace()?;
    let expected = GridFn::from_range_values((-2.0, 0.0), array![3.0, 2.0, 1.0])?;
    assert_eq!(expected, grid_fn);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_to_grid_finer() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![0.0, 10.0, 20.0])?;
    let resampled = grid_fn.resample_range_dx((0.0, 2.0), 0.5)?;
    assert_ulps_eq!(resampled.x_min(), 0.0);
    assert_ulps_eq!(resampled.x_max(), 2.0);
    assert_ulps_eq!(resampled.dx(), 0.5);
    assert_eq!(resampled.n_points(), 5);
    assert_ulps_eq!(resampled.y(), &array![0.0, 5.0, 10.0, 15.0, 20.0], max_ulps = 4);
    Ok(())
  }

  /// A soft boundary continues the distribution past the grid edge (`Linear`); a hard
  /// boundary terminates the domain (`Hard`: zero beyond; `Error`: undefined beyond).
  #[rustfmt::skip]
  #[rstest]
  #[case::linear(BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 }), true)]
  #[case::hard(  BoundaryBehavior::Hard,                                false)]
  #[case::error( BoundaryBehavior::Error,                               false)]
  #[trace]
  fn test_boundary_behavior_is_soft(#[case] behavior: BoundaryBehavior, #[case] expected: bool) {
    assert_eq!(expected, behavior.is_soft());
  }
}
