#[cfg(test)]
mod tests {
  use crate::*;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use rstest::rstest;
  use treetime_utils::assert_error;

  /// Build a GridFn whose stored neg-log ordinates follow a straight line `y(t) = slope * t`, for
  /// testing that `SoftTailLaw::fit` recovers the known neg-log slope. The fit regresses the
  /// stored ordinates against `t` directly, so it recovers `slope` exactly (subject to the decay
  /// guard).
  fn make_neglog_linear_grid(slope: f64, x_min: f64, x_max: f64, n: usize) -> GridFn<f64> {
    let y = Array1::linspace(x_min, x_max, n).mapv(|t| slope * t);
    GridFn::from_range_values((x_min, x_max), y).unwrap()
  }

  // --- SoftTailLaw::eval ---
  // Oracle: neg-log line `y(t) = y_edge + slope * (t - t_edge)`, evaluated by hand.

  #[rustfmt::skip]
  #[rstest]
  // Continuity: at the edge the tail equals the edge ordinate for any slope.
  #[case::at_edge_flat(    0.0, 2.0, 3.0, 3.0,  2.0)]
  #[case::at_edge_decaying(0.8, 2.0, 3.0, 3.0,  2.0)]
  // Right tail decays away from support (slope > 0): y_edge + slope * (t - t_edge).
  #[case::right_decay_unit(1.0, 2.0, 0.0, 1.0,  3.0)]
  #[case::right_decay_two( 0.5, 4.0, 0.0, 2.0,  5.0)]
  // Left tail decays away from support (slope < 0): outward is decreasing t.
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
    assert_abs_diff_eq!(expected, law.eval(GridEdge { t: t_edge, y: y_edge }, t), epsilon = 1e-14);
  }

  // --- SoftTailLaw::fit ---

  #[rustfmt::skip]
  #[rstest]
  // Right tail: stored neg-log line rises rightward (slope > 0), so p decays rightward. Recovered.
  #[case::right_decay_slow( 0.8, Side::Right,  0.8)]
  #[case::right_decay_fast( 2.5, Side::Right,  2.5)]
  // Left tail: stored neg-log line falls rightward (slope < 0), so p decays leftward. Recovered.
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
    assert_abs_diff_eq!(expected_slope, law.slope, epsilon = 1e-10);
  }

  #[test]
  fn test_soft_tail_law_fit_clamps_growing_right_tail_to_flat() {
    // A neg-log line that falls rightward (slope < 0) means p grows toward the right edge, so the
    // density does not decay away from support. The right-side decay guard clamps the slope to zero.
    let grid = make_neglog_linear_grid(-1.0, 0.1, 2.0, 20);
    let law = SoftTailLaw::fit(&grid, Side::Right, 10).expect("fit should succeed");
    assert_abs_diff_eq!(0.0, law.slope, epsilon = 1e-14);
  }

  #[test]
  fn test_soft_tail_law_fit_clamps_growing_left_tail_to_flat() {
    // A neg-log line that rises rightward (slope > 0) means p grows toward the left edge; the
    // left-side decay guard clamps the slope to zero.
    let grid = make_neglog_linear_grid(1.0, 0.1, 2.0, 20);
    let law = SoftTailLaw::fit(&grid, Side::Left, 10).expect("fit should succeed");
    assert_abs_diff_eq!(0.0, law.slope, epsilon = 1e-14);
  }

  #[test]
  fn test_soft_tail_law_fit_err_for_non_finite() -> Result<(), Report> {
    // Zero-probability points store as +inf under NegLog; a fit needs two finite ordinates.
    let grid = GridFn::from_range_values((0.1, 1.0), Array1::from_elem(10, f64::INFINITY))?;
    let law = SoftTailLaw::fit(&grid, Side::Right, 5);
    assert_error!(
      law,
      "Soft-tail fit on the Right side needs at least two finite grid points near the edge, found 0"
    );
    Ok(())
  }

  // --- SoftTailLaw::mass ---
  // Oracle: half-line integral exp(-y_edge) / |slope|, with p_edge = exp(-y_edge).

  #[rustfmt::skip]
  #[rstest]
  #[case::unit_slope_unit_edge( 1.0,  0.0,                    1.0)]
  #[case::half_slope(           0.5,  0.0,                    2.0)]
  #[case::negative_slope_abs(  -0.5,  0.0,                    2.0)]
  #[case::nonzero_edge(         1.0,  std::f64::consts::LN_2, 0.5)]
  #[trace]
  fn test_soft_tail_law_mass(#[case] slope: f64, #[case] y_edge: f64, #[case] expected: f64) {
    let law = SoftTailLaw { slope };
    assert_abs_diff_eq!(expected, law.mass(y_edge), epsilon = 1e-14);
  }

  #[test]
  fn test_soft_tail_law_mass_flat_is_infinite() {
    // A flat tail (slope 0) is non-integrable, unlike a genuine decaying tail.
    let law = SoftTailLaw { slope: 0.0 };
    assert!(law.mass(2.0).is_infinite());
  }

  // --- SoftTailLaw::compose_multiply ---

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
    // Multiplication is addition in neg-log space, so the composed slope is the sum.
    let composed = SoftTailLaw { slope: slope_a }.compose_multiply(&SoftTailLaw { slope: slope_b });
    assert_abs_diff_eq!(expected, composed.slope, epsilon = 1e-14);
  }

  // --- SoftTailLaw::negate_arg ---

  #[rstest]
  #[case::positive(0.8)]
  #[case::negative(-1.3)]
  #[case::zero(0.0)]
  #[trace]
  fn test_soft_tail_law_negate_arg_flips_sign_and_is_involution(#[case] slope: f64) {
    let law = SoftTailLaw { slope };
    assert_abs_diff_eq!(-slope, law.negate_arg().slope, epsilon = 1e-14);
    // Negating twice restores the original slope.
    assert_abs_diff_eq!(slope, law.negate_arg().negate_arg().slope, epsilon = 1e-14);
  }

  // --- GridFn integration: soft tail in extrapolation ---

  #[test]
  fn test_gridfn_soft_tail_right_extrapolation() -> Result<(), Report> {
    // Grid [0.0, 2.0], soft right tail with neg-log slope 1.0 (decaying rightward).
    let grid = GridFn::from_range_values((0.0, 2.0), ndarray::array![4.0, 3.0, 2.0])?;
    let right = BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 });

    // At the edge: continuous with the grid.
    assert_abs_diff_eq!(2.0, grid.interp(2.0)?, epsilon = 1e-14);
    // Beyond the edge: y_edge + slope * (t - x_max) = 2.0 + 1.0 * 0.5.
    assert_abs_diff_eq!(
      2.5,
      grid.interp_with_extrap(2.5, BoundaryBehavior::Error, right)?,
      epsilon = 1e-14
    );
    Ok(())
  }

  #[test]
  fn test_gridfn_soft_tail_left_extrapolation() -> Result<(), Report> {
    // Grid [0.0, 2.0], soft left tail with neg-log slope -1.0 (decaying leftward).
    let grid = GridFn::from_range_values((0.0, 2.0), ndarray::array![2.0, 3.0, 4.0])?;
    let left = BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 });

    assert_abs_diff_eq!(2.0, grid.interp(0.0)?, epsilon = 1e-14);
    // Beyond the left edge: y_edge + slope * (t - x_min) = 2.0 + (-1.0) * (-0.5).
    assert_abs_diff_eq!(
      2.5,
      grid.interp_with_extrap(-0.5, left, BoundaryBehavior::Error)?,
      epsilon = 1e-14
    );
    Ok(())
  }

  #[test]
  fn test_gridfn_soft_tail_is_edge_relative_after_rewindowing() -> Result<(), Report> {
    // A soft edge is movable: re-windowing to a narrower grid moves the edge inward. Because the
    // tail is anchored to the live edge (not an absolute coordinate), extrapolation past the new
    // edge uses the new edge value, staying continuous. An absolutely anchored law would extend
    // from the stale original edge and be discontinuous here.
    let grid = make_neglog_linear_grid(-1.0, 0.0, 5.0, 51);
    let right = BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 });

    // Re-window to [0, 4]: the right edge moves from 5.0 to 4.0.
    let rewindowed = grid.resample_range_dx((0.0, 4.0), 0.1)?;
    let edge_value = rewindowed.interp(4.0)?;
    // Extrapolation just past the new edge equals edge_value + slope * 0.5.
    assert_abs_diff_eq!(
      edge_value + 0.5,
      rewindowed.interp_with_extrap(4.5, BoundaryBehavior::Error, right)?,
      epsilon = 1e-12
    );
    Ok(())
  }
}
