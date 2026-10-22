#[cfg(test)]
mod tests {
  use crate::*;
  use eyre::Report;
  use ndarray::Array1;
  use rstest::rstest;
  use treetime_utils::pretty_assert_abs_diff_eq;
  use treetime_utils::{assert_error, pretty_assert_ulps_eq};

  fn make_neglog_linear_grid(slope: f64, x_min: f64, x_max: f64, n: usize) -> GridFn<f64> {
    let y = Array1::linspace(x_min, x_max, n).mapv(|t| slope * t);
    GridFn::from_range_values((x_min, x_max), y).unwrap()
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::at_edge_flat(    0.0, 2.0, 3.0, 3.0,  2.0)]
  #[case::at_edge_decaying(0.8, 2.0, 3.0, 3.0,  2.0)]
  #[case::right_decay_unit(1.0, 2.0, 0.0, 1.0,  3.0)]
  #[case::right_decay_two( 0.5, 4.0, 0.0, 2.0,  5.0)]
  #[case::left_decay_unit(-1.0, 2.0, 0.0, -1.0, 3.0)]
  #[trace]
  fn test_soft_tail_law_eval(
    #[case] slope: f64,
    #[case] y_edge: f64,
    #[case] t_edge: f64,
    #[case] t: f64,
    #[case] expected: f64,
  ) {
    let law = SoftTailLaw { slope };
    pretty_assert_abs_diff_eq!(expected, law.eval(GridEdge { t: t_edge, y: y_edge }, t), epsilon = 1e-14);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::right_decay_slow( 0.8, Side::Right,  0.8)]
  #[case::right_decay_fast( 2.5, Side::Right,  2.5)]
  #[case::left_decay_slow( -0.6, Side::Left,  -0.6)]
  #[case::left_decay_fast( -3.1, Side::Left,  -3.1)]
  #[trace]
  fn test_soft_tail_law_fit_recovers_exact_slope(
    #[case] slope: f64,
    #[case] side: Side,
    #[case] expected_slope: f64,
  ) {
    let grid = make_neglog_linear_grid(slope, 0.1, 2.0, 20);
    let law = SoftTailLaw::fit(&grid, side, 10).expect("fit should succeed");
    pretty_assert_abs_diff_eq!(expected_slope, law.slope, epsilon = 1e-10);
  }

  #[test]
  fn test_soft_tail_law_fit_clamps_growing_right_tail_to_flat() {
    let grid = make_neglog_linear_grid(-1.0, 0.1, 2.0, 20);
    let law = SoftTailLaw::fit(&grid, Side::Right, 10).expect("fit should succeed");
    pretty_assert_abs_diff_eq!(0.0, law.slope, epsilon = 1e-14);
  }

  #[test]
  fn test_soft_tail_law_fit_clamps_growing_left_tail_to_flat() {
    let grid = make_neglog_linear_grid(1.0, 0.1, 2.0, 20);
    let law = SoftTailLaw::fit(&grid, Side::Left, 10).expect("fit should succeed");
    pretty_assert_abs_diff_eq!(0.0, law.slope, epsilon = 1e-14);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::right_decay(       ( 2.0, Side::Right), (1.0, 0.5), 1.5,  1.5)]
  #[case::left_decay(        (-0.5, Side::Left),  (1.0, 0.5), 1.5, -1.0)]
  #[case::right_far_level(   ( 1e-3, Side::Right),(0.0, 0.0), 30.0, 3e4)]
  #[case::at_edge_ordinate(  ( 2.0, Side::Right), (1.0, 0.5), 0.5,  1.0)]
  #[case::below_edge_ordinate((-2.0, Side::Left), (1.0, 0.5), 0.1,  1.0)]
  #[trace]
  fn test_soft_tail_law_position_at_inverts_eval(
    #[case] (slope, side): (f64, Side),
    #[case] (t_edge, y_edge): (f64, f64),
    #[case] level: f64,
    #[case] expected: f64,
  ) {
    let law = SoftTailLaw { slope };
    let edge = GridEdge { t: t_edge, y: y_edge };
    let position = law.position_at(edge, level, side);
    pretty_assert_abs_diff_eq!(expected, position, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(level.max(y_edge), law.eval(edge, position), epsilon = 1e-15);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::flat_right(                ( 0.0, Side::Right), 2.0,  f64::INFINITY)]
  #[case::flat_left(                 ( 0.0, Side::Left),  2.0, -f64::INFINITY)]
  #[case::flat_at_edge_ordinate(     ( 0.0, Side::Right), 0.5,  f64::INFINITY)]
  #[case::flat_below_edge_ordinate(  ( 0.0, Side::Right), 0.4,  1.0)]
  #[case::growing_right(             (-1.0, Side::Right), 2.0,  f64::INFINITY)]
  #[case::growing_left(              ( 1.0, Side::Left),  2.0, -f64::INFINITY)]
  #[trace]
  fn test_soft_tail_law_position_at_level_held_without_decay_is_unbounded(
    #[case] (slope, side): (f64, Side),
    #[case] level: f64,
    #[case] expected: f64,
  ) {
    let law = SoftTailLaw { slope };
    pretty_assert_ulps_eq!(expected, law.position_at(GridEdge { t: 1.0, y: 0.5 }, level, side), max_ulps = 0);
  }

  #[test]
  fn test_soft_tail_law_fit_err_for_non_finite() -> Result<(), Report> {
    let grid = GridFn::from_range_values((0.1, 1.0), Array1::from_elem(10, f64::INFINITY))?;
    let law = SoftTailLaw::fit(&grid, Side::Right, 5);
    assert_error!(
      law,
      "Soft-tail fit on the Right side needs at least two finite grid points near the edge, found 0"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::unit_slope_unit_edge( 1.0,  0.0,                    1.0)]
  #[case::half_slope(           0.5,  0.0,                    2.0)]
  #[case::negative_slope_abs(  -0.5,  0.0,                    2.0)]
  #[case::nonzero_edge(         1.0,  std::f64::consts::LN_2, 0.5)]
  #[trace]
  fn test_soft_tail_law_mass(#[case] slope: f64, #[case] y_edge: f64, #[case] expected: f64) {
    let law = SoftTailLaw { slope };
    pretty_assert_abs_diff_eq!(expected, law.mass(y_edge), epsilon = 1e-14);
  }

  #[test]
  fn test_soft_tail_law_mass_flat_is_infinite() {
    let law = SoftTailLaw { slope: 0.0 };
    assert!(law.mass(2.0).is_infinite());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::both_decay(     0.8,  1.5,  2.3)]
  #[case::opposite_signs( 1.0, -0.4,  0.6)]
  #[case::with_flat(      0.7,  0.0,  0.7)]
  #[trace]
  fn test_soft_tail_law_compose_multiply_slopes_add(
    #[case] slope_a: f64,
    #[case] slope_b: f64,
    #[case] expected: f64,
  ) {
    let composed = SoftTailLaw { slope: slope_a }.compose_multiply(&SoftTailLaw { slope: slope_b });
    pretty_assert_abs_diff_eq!(expected, composed.slope, epsilon = 1e-14);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::reflect_positive( 0.8, -1.0, -0.8)]
  #[case::reflect_negative(-1.3, -1.0,  1.3)]
  #[case::reflect_flat(     0.0, -1.0,  0.0)]
  #[case::stretch(          1.5,  2.0,  0.75)]
  #[case::compress(        -1.5,  0.5, -3.0)]
  #[trace]
  fn test_soft_tail_law_scale_arg_divides_the_slope(#[case] slope: f64, #[case] factor: f64, #[case] expected: f64) {
    pretty_assert_abs_diff_eq!(expected, SoftTailLaw { slope }.scale_arg(factor).slope, epsilon = 1e-14);
  }

  #[test]
  fn test_soft_tail_law_reflection_is_involution() {
    let law = SoftTailLaw { slope: 0.8 };
    assert_eq!(law, law.scale_arg(-1.0).scale_arg(-1.0));
  }

  #[test]
  fn test_gridfn_soft_tail_right_extrapolation() -> Result<(), Report> {
    let grid = GridFn::from_range_values((0.0, 2.0), ndarray::array![4.0, 3.0, 2.0])?;
    let right = BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 });

    pretty_assert_abs_diff_eq!(2.0, grid.interp(2.0)?, epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(
      2.5,
      grid.interp_with_extrap(2.5, BoundaryBehavior::Error, right)?,
      epsilon = 1e-14
    );
    Ok(())
  }

  #[test]
  fn test_gridfn_soft_tail_left_extrapolation() -> Result<(), Report> {
    let grid = GridFn::from_range_values((0.0, 2.0), ndarray::array![2.0, 3.0, 4.0])?;
    let left = BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 });

    pretty_assert_abs_diff_eq!(2.0, grid.interp(0.0)?, epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(
      2.5,
      grid.interp_with_extrap(-0.5, left, BoundaryBehavior::Error)?,
      epsilon = 1e-14
    );
    Ok(())
  }
}
