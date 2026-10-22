#[cfg(test)]
mod tests {
  use crate::MaxGridPoints;
  use crate::*;
  use eyre::Report;
  use ndarray::{Array1, array, s};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::ops::RangeInclusive;
  use treetime_utils::assert_error;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

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
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 4);
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
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_range_dx_clamped_holds_boundary_on_overshoot() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, 10.0])?;
    let clamped = grid_fn.resample_range_dx_clamped((0.0, 1.0), 0.4, MaxGridPoints::default())?;
    pretty_assert_ulps_eq!(clamped.y(), &array![0.0, 4.0, 8.0, 10.0], max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_range_dx_clamped_interpolates_within_support() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, 10.0])?;
    let clamped = grid_fn.resample_range_dx_clamped((0.0, 1.0), 0.5, MaxGridPoints::default())?;
    pretty_assert_ulps_eq!(clamped.y(), &array![0.0, 5.0, 10.0], max_ulps = 4);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
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
    pretty_assert_ulps_eq!(0.0, actual, max_ulps = 4);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::below_a_unit_grid_within_rounding(   (0.0,    1.0),    -4e-16,         2.0)]
  #[case::below_a_unit_grid_beyond_rounding(   (0.0,    1.0),    -2e-15,         0.0)]
  #[case::below_a_small_grid_beyond_rounding(  (0.0,    1e-6),   -1e-17,         0.0)]
  #[case::above_a_date_grid_within_rounding(   (2000.0, 2010.0), 2010.0 + 1e-12, 10.0)]
  #[case::above_a_date_grid_beyond_rounding(   (2000.0, 2010.0), 2010.0 + 1e-11, 0.0)]
  #[trace]
  fn test_gridfn_interp_just_outside_the_grid_is_the_boundary_value_within_rounding_of_the_grid_magnitude(
    #[case] x_range: (f64, f64),
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values(x_range, array![2.0, 10.0])?;
    let actual = grid_fn.interp_with_extrap(query, BoundaryBehavior::Hard, BoundaryBehavior::Hard)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::left(-1.0, "GridFn evaluated at -1.0, below the support boundary 0.0, but no extrapolation policy is set for that side")]
  #[case::right(2.0, "GridFn evaluated at 2.0, above the support boundary 1.0, but no extrapolation policy is set for that side")]
  #[trace]
  fn test_gridfn_interp_out_of_support_errors_by_default(#[case] query: f64, #[case] expected: &str) -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![2.0, 10.0])?;
    assert_error!(grid_fn.interp(query), expected);
    Ok(())
  }

  #[test]
  fn test_gridfn_accessors() -> Result<(), Report> {
    let x_range = (0.0, 2.0);
    let y = array![10.0, 20.0, 30.0];
    let grid_fn = GridFn::from_range_values(x_range, y.clone())?;
    let x = grid_fn.x();
    pretty_assert_ulps_eq!(x[0], 0.0);
    pretty_assert_ulps_eq!(x[1], 1.0);
    pretty_assert_ulps_eq!(x[2], 2.0);
    assert_eq!(&y, grid_fn.y());
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

  #[rustfmt::skip]
  #[rstest]
  #[case::linear(BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 }), true)]
  #[case::hard(  BoundaryBehavior::Hard,                                false)]
  #[case::error( BoundaryBehavior::Error,                               false)]
  #[trace]
  fn test_boundary_behavior_is_soft(#[case] behavior: BoundaryBehavior, #[case] expected: bool) {
    assert_eq!(expected, behavior.is_soft());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::cell_midpoint(  9.5)]
  #[case::near_boundary(  9.99)]
  #[case::near_anchor(    9.01)]
  #[trace]
  fn test_gridfn_interp_boundary_cell_follows_the_hard_approach_law(#[case] query: f64) -> Result<(), Report> {
    let law = |t: f64| 1.0 - 2.0 * (10.0 - t).ln();
    let t = Array1::linspace(0.0, 10.0, 11);
    let y = t.mapv(|ti: f64| if ti < 10.0 { law(ti) } else { f64::INFINITY });
    let grid_fn = GridFn::from_range_values((0.0, 10.0), y)?;

    let actual = grid_fn.interp(query)?;

    pretty_assert_abs_diff_eq!(law(query), actual, epsilon = 1e-9);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_boundary_cell_uses_the_stored_law_of_its_side() -> Result<(), Report> {
    let y = array![0.0, 0.0, 0.0, 0.0, f64::INFINITY];
    let grid_fn = GridFn::from_range_values((0.0, 4.0), y)?;
    let law = BoundaryBehavior::HardApproach(HardApproachLaw::from_points(
      4.0,
      Side::Right,
      &[(3.0, 0.0), (2.0, -3.0 * 2.0_f64.ln())],
    )?);

    let actual = grid_fn.interp_with_extrap(3.5, BoundaryBehavior::Hard, law)?;

    pretty_assert_abs_diff_eq!(-3.0 * 0.5_f64.ln(), actual, epsilon = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::grid_point_next_to_the_boundary(  2.0 + 1e-13, 0.0)]
  #[case::left_end_infinite(                0.5,         f64::NAN)]
  #[case::both_ends_infinite(               3.5,         f64::INFINITY)]
  #[trace]
  fn test_gridfn_interp_never_returns_nan_next_to_zero_probability(
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let y = array![f64::INFINITY, 2.0, 0.0, f64::INFINITY, f64::INFINITY];
    let grid_fn = GridFn::from_range_values((0.0, 4.0), y)?;

    let actual = grid_fn.interp(query)?;

    assert!(!actual.is_nan(), "interpolation at {query} returned NaN");
    if expected.is_infinite() {
      assert!(actual.is_infinite() && actual > 0.0, "expected +inf at {query}, got {actual}");
    } else if !expected.is_nan() {
      pretty_assert_ulps_eq!(expected, actual, max_ulps = 4);
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::regular_cell(        array![0.0, 1.0, 2.0],                     0, BoundaryBehavior::Hard, None)]
  #[case::both_ends_infinite(  array![0.0, f64::INFINITY, f64::INFINITY], 1, BoundaryBehavior::Hard, None)]
  #[case::flat_without_points( array![f64::INFINITY, 2.0, f64::INFINITY], 0, BoundaryBehavior::Hard,
                               Some(HardCell { law: HardApproachLaw { t_hard: 0.0, b: 0.0 }, anchor: GridEdge { t: 1.0, y: 2.0 } }))]
  #[case::stored_law_of_side(  array![0.0, 0.0, 0.0, 0.0, f64::INFINITY], 3, BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 4.0, b: 1.5 }),
                               Some(HardCell { law: HardApproachLaw { t_hard: 4.0, b: 1.5 }, anchor: GridEdge { t: 3.0, y: 0.0 } }))]
  #[case::stored_law_elsewhere(array![0.0, 0.0, 0.0, 0.0, f64::INFINITY], 3, BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 5.0, b: 1.5 }),
                               Some(HardCell { law: HardApproachLaw { t_hard: 4.0, b: 0.0 }, anchor: GridEdge { t: 3.0, y: 0.0 } }))]
  #[trace]
  fn test_gridfn_hard_cell_returns_the_law_of_a_cell_with_one_infinite_end(
    #[case] y: Array1<f64>,
    #[case] cell: usize,
    #[case] right: BoundaryBehavior,
    #[case] expected: Option<HardCell>,
  ) -> Result<(), Report> {
    let n = y.len().to_f64_exact();
    let grid_fn = GridFn::from_range_values((0.0, n - 1.0), y)?;
    let actual = grid_fn.hard_cell(cell, BoundaryBehavior::Hard, right)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_hard_cell_fits_the_power_law_of_the_finite_points() -> Result<(), Report> {
    let y = array![
      -3.0 * 4.0_f64.ln(),
      -3.0 * 3.0_f64.ln(),
      -3.0 * 2.0_f64.ln(),
      0.0,
      f64::INFINITY
    ];
    let grid_fn = GridFn::from_range_values((0.0, 4.0), y)?;
    let actual = grid_fn
      .hard_cell(3, BoundaryBehavior::Hard, BoundaryBehavior::Hard)?
      .unwrap();
    assert_eq!(GridEdge { t: 3.0, y: 0.0 }, actual.anchor);
    pretty_assert_ulps_eq!(4.0, actual.law.t_hard, max_ulps = 0);
    pretty_assert_ulps_eq!(3.0, actual.law.b, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_gridfn_hard_cell_value_is_the_law_of_the_cell() -> Result<(), Report> {
    let y = array![f64::INFINITY, 1.0, 1.5, 2.5, 4.0, f64::INFINITY];
    let grid_fn = GridFn::from_range_values((0.0, 5.0), y)?;
    for (cell, t) in [(0, 0.3), (4, 4.6)] {
      let hard_cell = grid_fn
        .hard_cell(cell, BoundaryBehavior::Hard, BoundaryBehavior::Hard)?
        .unwrap();
      pretty_assert_ulps_eq!(hard_cell.eval(t), grid_fn.interp(t)?, max_ulps = 0);
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::right_end_reads_five_points(  array![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, f64::INFINITY], 5, Some(1..=5))]
  #[case::right_end_stops_at_infinity(  array![0.0, f64::INFINITY, 0.0, 0.0, f64::INFINITY], 3, Some(2..=3))]
  #[case::right_end_stops_at_grid_start(array![0.0, 0.0, f64::INFINITY],                     1, Some(0..=1))]
  #[case::left_end_reads_five_points(   array![f64::INFINITY, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 0, Some(1..=5))]
  #[case::left_end_stops_at_grid_end(   array![1.0, f64::INFINITY, 2.0, 3.0],                1, Some(2..=3))]
  #[case::finite_cell(                  array![0.0, 1.0],                                    0, None)]
  #[case::zero_cell(                    array![f64::INFINITY, f64::INFINITY],                0, None)]
  #[trace]
  fn test_gridfn_hard_cell_fit_range_is_the_points_the_fitted_law_reads(
    #[case] y: Array1<f64>,
    #[case] cell: usize,
    #[case] expected: Option<RangeInclusive<usize>>,
  ) -> Result<(), Report> {
    let n = y.len().to_f64_exact();
    let grid_fn = GridFn::from_range_values((0.0, n - 1.0), y)?;
    assert_eq!(expected, grid_fn.hard_cell_fit_range(cell)?);
    Ok(())
  }

  #[test]
  fn test_gridfn_hard_cell_law_depends_only_on_its_fit_range() -> Result<(), Report> {
    let y = array![7.0, 30.0, 30.0, 0.0, 300.0, 0.0, f64::INFINITY];
    let grid_fn = GridFn::from_range_values((0.0, 6.0), y.clone())?;
    let range = grid_fn.hard_cell_fit_range(5)?.unwrap();
    let sliced = GridFn::from_range_values((1.0, 6.0), y.slice(s![*range.start()..]).to_owned())?;
    let full = grid_fn.hard_cell(5, BoundaryBehavior::Hard, BoundaryBehavior::Hard)?;
    let kept = sliced.hard_cell(4, BoundaryBehavior::Hard, BoundaryBehavior::Hard)?;
    assert_eq!(1..=5, range);
    assert_eq!(full, kept);
    Ok(())
  }

  #[test]
  fn test_gridfn_hard_cell_fit_range_outside_the_grid_is_error() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, f64::INFINITY])?;
    assert_error!(grid_fn.hard_cell_fit_range(1), "Cell 1 is outside the grid of 2 points");
    Ok(())
  }

  #[test]
  fn test_gridfn_hard_cell_outside_the_grid_is_error() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, f64::INFINITY])?;
    assert_error!(
      grid_fn.hard_cell(1, BoundaryBehavior::Hard, BoundaryBehavior::Hard),
      "Cell 1 is outside the grid of 2 points"
    );
    Ok(())
  }
}
