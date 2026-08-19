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
    let grid_fn = GridFn::from_range_values(x_range, y)?.with_extrap(BoundaryBehavior::Hard);
    let actual = grid_fn.interp(query)?;
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
  fn test_gridfn_negate_arg_swaps_extrap_sides() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?
      .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 }))
      .with_right_extrap(BoundaryBehavior::Hard);
    let negated = grid_fn.negate_arg()?;
    // Sides swap; the reflected soft law's slope flips sign on its new (right) side.
    assert_eq!(BoundaryBehavior::Hard, negated.left_extrap());
    assert_eq!(
      BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 }),
      negated.right_extrap()
    );
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_preserves_extrap() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![0.0, 10.0, 20.0])?
      .with_left_extrap(BoundaryBehavior::Hard)
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 2.0 }));
    let resampled = grid_fn.resample_range_dx((0.0, 2.0), 0.5)?;
    assert_eq!(BoundaryBehavior::Hard, resampled.left_extrap());
    assert_eq!(
      BoundaryBehavior::Linear(SoftTailLaw { slope: 2.0 }),
      resampled.right_extrap()
    );
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

  /// A raw grid and y-array carry no declared out-of-support policy, so both tails must default to
  /// `Error` on fresh construction. Regridding an existing function is the case that must instead
  /// preserve policy; this test pins the deliberate fresh-construction reset so the two paths stay
  /// distinct.
  #[test]
  fn test_gridfn_from_grid_array_defaults_to_error() -> Result<(), Report> {
    let grid = Grid::from_range_n_points(0.0, 2.0, 3)?;
    let grid_fn = GridFn::from_grid_array(grid, array![0.0, 1.0, 2.0])?;
    assert_eq!(BoundaryBehavior::Error, grid_fn.left_extrap());
    assert_eq!(BoundaryBehavior::Error, grid_fn.right_extrap());
    Ok(())
  }

  /// Regridding must carry the fitted `HardApproachLaw`, not only the hard/soft class. The law is
  /// edge-relative (stores only the exponent), so a pure resample preserves it verbatim.
  #[test]
  fn test_gridfn_resample_preserves_fitted_hard_law() -> Result<(), Report> {
    let law = HardApproachLaw {
      t_hard: 0.0,
      shape: Approach::Divergent { b: 1.0 },
    };
    let grid_fn = GridFn::from_range_values((1.0, 3.0), array![2.0, 3.0, 4.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(law));

    let resampled = grid_fn.resample_range_dx((1.0, 3.0), 0.25)?;

    assert_eq!(BoundaryBehavior::HardApproach(law), resampled.left_extrap());
    assert!(
      !resampled.left_extrap().is_soft(),
      "a hard boundary must stay hard across regridding"
    );
    Ok(())
  }

  /// Regridding must carry the fitted `SoftTailLaw`. The soft-tail law is edge-relative and reads
  /// the resampled edge ordinate on evaluation, so a pure resample preserves its slope verbatim.
  #[test]
  fn test_gridfn_resample_preserves_fitted_soft_law() -> Result<(), Report> {
    let law = SoftTailLaw { slope: 0.7 };
    let grid_fn =
      GridFn::from_range_values((1.0, 3.0), array![4.0, 3.0, 2.0])?.with_right_extrap(BoundaryBehavior::Linear(law));

    let resampled = grid_fn.resample_range_dx((1.0, 3.0), 0.25)?;

    assert_eq!(BoundaryBehavior::Linear(law), resampled.right_extrap());
    assert!(
      resampled.right_extrap().is_soft(),
      "a soft boundary must stay soft across regridding"
    );
    Ok(())
  }

  /// Scaling every neg-log ordinate by a factor raises the probability to a power (`p -> p^factor`),
  /// so the hard approach law's shape parameter (the exponent `b`) scales by the same factor while
  /// the boundary location is unchanged. The hard class must survive the scale.
  #[test]
  fn test_gridfn_scale_y_preserves_hard_class_and_scales_law() -> Result<(), Report> {
    let law = HardApproachLaw {
      t_hard: 0.0,
      shape: Approach::Divergent { b: 1.0 },
    };
    let grid_fn = GridFn::from_range_values((1.0, 3.0), array![2.0, 3.0, 4.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(law));

    let scaled = grid_fn.scale_y(3.0);

    let expected = HardApproachLaw {
      t_hard: 0.0,
      shape: Approach::Divergent { b: 3.0 },
    };
    assert_eq!(BoundaryBehavior::HardApproach(expected), scaled.left_extrap());
    Ok(())
  }

  /// Scaling every neg-log ordinate by a factor scales the soft-tail slope by the same factor. The
  /// soft class must survive the scale.
  #[test]
  fn test_gridfn_scale_y_preserves_soft_class_and_scales_law() -> Result<(), Report> {
    let law = SoftTailLaw { slope: 0.7 };
    let grid_fn =
      GridFn::from_range_values((1.0, 3.0), array![4.0, 3.0, 2.0])?.with_right_extrap(BoundaryBehavior::Linear(law));

    let scaled = grid_fn.scale_y(2.0);

    let expected = SoftTailLaw { slope: 0.7 * 2.0 };
    assert_eq!(BoundaryBehavior::Linear(expected), scaled.right_extrap());
    Ok(())
  }

  /// Adding a constant to every neg-log ordinate (a `NegLog` normalize shift) leaves the hard
  /// approach law unchanged: it is edge-relative and stores only the shift-invariant exponent, so
  /// it reads the shifted edge ordinate on evaluation. The hard class must survive the shift.
  #[test]
  fn test_gridfn_shift_y_preserves_hard_law_unchanged() -> Result<(), Report> {
    let law = HardApproachLaw {
      t_hard: 0.0,
      shape: Approach::Divergent { b: 1.0 },
    };
    let grid_fn = GridFn::from_range_values((1.0, 3.0), array![2.0, 3.0, 4.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(law));

    let shifted = grid_fn.shift_y(-3.0);

    assert_eq!(BoundaryBehavior::HardApproach(law), shifted.left_extrap());
    Ok(())
  }

  /// Adding a constant to every neg-log ordinate leaves a soft-tail slope unchanged: the slope is
  /// shift-invariant and the edge-relative tail reads the shifted edge ordinate on evaluation. The
  /// soft class and the fitted law must survive the shift verbatim.
  #[test]
  fn test_gridfn_shift_y_preserves_soft_law_unchanged() -> Result<(), Report> {
    let law = SoftTailLaw { slope: 0.7 };
    let grid_fn =
      GridFn::from_range_values((1.0, 3.0), array![4.0, 3.0, 2.0])?.with_right_extrap(BoundaryBehavior::Linear(law));

    let shifted = grid_fn.shift_y(1.5);

    assert_eq!(BoundaryBehavior::Linear(law), shifted.right_extrap());
    Ok(())
  }

  /// End-to-end round trip through `from_grid_array`-based regridding (`resample`) and `scale_y`:
  /// a hard boundary stays hard, a soft boundary stays soft, and both fitted laws survive intact.
  /// `scale_y(1.0)` is the identity, so every coefficient is preserved exactly across the trip.
  #[test]
  fn test_gridfn_boundary_law_survives_regrid_round_trip() -> Result<(), Report> {
    let hard = HardApproachLaw {
      t_hard: 0.0,
      shape: Approach::Divergent { b: 1.0 },
    };
    let soft = SoftTailLaw { slope: 0.7 };
    let grid_fn = GridFn::from_range_values((1.0, 3.0), array![3.0, 2.5, 2.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(hard))
      .with_right_extrap(BoundaryBehavior::Linear(soft));

    let round_tripped = grid_fn
      .resample_range_dx((1.0, 3.0), 0.25)?
      .scale_y(1.0)
      .resample_range_dx((1.0, 3.0), 0.5)?;

    assert_eq!(BoundaryBehavior::HardApproach(hard), round_tripped.left_extrap());
    assert_eq!(BoundaryBehavior::Linear(soft), round_tripped.right_extrap());
    assert!(!round_tripped.left_extrap().is_soft(), "hard boundary must stay hard");
    assert!(round_tripped.right_extrap().is_soft(), "soft boundary must stay soft");
    Ok(())
  }
}
