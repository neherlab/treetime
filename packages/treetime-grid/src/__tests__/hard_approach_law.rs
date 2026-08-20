#[cfg(test)]
mod tests {
  use crate::*;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  /// Build a GridFn whose stored ordinates follow neg-log of a power law:
  /// `y(t) = a - b * ln|t - t_hard|`, for testing that the fit recovers the exponent `b`.
  fn make_neglog_power_law_grid(t_hard: f64, a: f64, b: f64, x_min: f64, x_max: f64, n: usize) -> GridFn<f64> {
    let y = Array1::linspace(x_min, x_max, n).mapv(|t| a - b * (t - t_hard).abs().ln());
    GridFn::from_range_values((x_min, x_max), y).unwrap()
  }

  // --- HardApproachLaw::eval (edge-relative) ---
  // Oracle: the closed form y(t) = edge.y - b * ln(|t - t_hard| / |edge.t - t_hard|), by hand.
  //   b > 0: diverges to +inf at t_hard (density zero there).
  //   b = 0: flat at edge.y everywhere including the boundary (finite density).

  #[rustfmt::skip]
  #[rstest]
  // Divergent b=1, edge (t_edge=1, y_edge=0): y(t) = -ln(t)
  #[case::b1_edge(       (0.0, 1.0), (0.0, 1.0), 1.0,                    0.0)]
  #[case::b1_half(       (0.0, 1.0), (0.0, 1.0), 0.5,          2.0_f64.ln())]
  #[case::b1_boundary(   (0.0, 1.0), (0.0, 1.0), 0.0,          f64::INFINITY)]
  // Divergent b=2, edge (1, 0): y(t) = -2 ln(t)
  #[case::b2_half(       (0.0, 2.0), (0.0, 1.0), 0.5,    2.0 * 2.0_f64.ln())]
  // Non-zero edge ordinate anchors the law: b=1, edge (t_edge=2, y_edge=3), y(1) = 3 - ln(1/2)
  #[case::anchored(      (0.0, 1.0), (3.0, 2.0), 1.0,      3.0 + 2.0_f64.ln())]
  // Flat b=0: value is the edge ordinate everywhere, including at the boundary
  #[case::flat_gap(      (0.0, 0.0), (2.5, 1.0), 0.3,                    2.5)]
  #[case::flat_boundary( (0.0, 0.0), (2.5, 1.0), 0.0,                    2.5)]
  // Non-zero t_hard, divergent b=1, edge (t_edge=6, y_edge=1): y(t) = 1 - ln(|t-5|/1)
  #[case::offset_edge(   (5.0, 1.0), (1.0, 6.0), 6.0,                    1.0)]
  #[case::offset_boundary((5.0, 1.0),(1.0, 6.0), 5.0,          f64::INFINITY)]
  #[case::offset_half(   (5.0, 1.0), (1.0, 6.0), 5.5,          1.0 + 2.0_f64.ln())]
  #[trace]
  fn test_hard_approach_law_eval(
    #[case] (t_hard, b): (f64, f64),
    #[case] (y_edge, t_edge): (f64, f64),
    #[case] t: f64,
    #[case] expected: f64,
  ) {
    let law = HardApproachLaw { t_hard, b };
    let actual = law.eval(GridEdge { t: t_edge, y: y_edge }, t);
    if expected.is_infinite() {
      assert!(actual.is_infinite() && actual > 0.0, "expected +inf, got {actual}");
    } else {
      assert_abs_diff_eq!(expected, actual, epsilon = 1e-14);
    }
  }

  // --- HardApproachLaw::fit (single-exponent power law) ---
  // Exact neg-log power-law data; the fit must recover the exponent `b` (the discarded intercept `a`
  // does not affect the edge-relative law).

  #[rustfmt::skip]
  #[rstest]
  #[case::b1(  0.0, 1.0, 1.0)]
  #[case::b2(  0.0, 0.5, 2.0)]
  #[case::b3(  0.0, 2.0, 3.0)]
  #[case::b05( 0.0, 1.0, 0.5)]
  #[trace]
  fn test_hard_approach_law_fit_recovers_exponent(
    #[case] t_hard: f64,
    #[case] a: f64,
    #[case] b: f64,
  ) -> Result<(), Report> {
    let grid = make_neglog_power_law_grid(t_hard, a, b, 0.1, 1.0, 20);
    let law = HardApproachLaw::fit(&grid, t_hard, Side::Left, 10).expect("fit should succeed");
    assert_abs_diff_eq!(t_hard, law.t_hard, epsilon = 1e-14);
    assert_abs_diff_eq!(b, law.b, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_hard_approach_law_fit_right_side() -> Result<(), Report> {
    // y(t) = 2.0 - 1.5 * ln(10 - t), divergent right boundary at t_hard = 10.0, fit from rightmost
    // points.
    let t_hard = 10.0;
    let b = 1.5;
    let y = Array1::linspace(8.0, 9.9, 20).mapv(|t: f64| 2.0 - b * (t_hard - t).ln());
    let grid = GridFn::from_range_values((8.0, 9.9), y)?;
    let law = HardApproachLaw::fit(&grid, t_hard, Side::Right, 10).expect("fit should succeed");
    assert_abs_diff_eq!(b, law.b, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_hard_approach_law_fit_clamps_wrong_sign_to_zero() -> Result<(), Report> {
    // A regression that trends the wrong way (density increasing into the boundary) clamps b to 0, a
    // flat shape, rather than a fabricated divergence. Ordinates decreasing with distance from t_hard
    // give a negative raw exponent, clamped to `b = 0`.
    let y = Array1::linspace(0.1, 1.0, 20).mapv(|t: f64| 2.0 + (t - 0.0).abs().ln());
    let grid = GridFn::from_range_values((0.1, 1.0), y)?;
    let law = HardApproachLaw::fit(&grid, 0.0, Side::Left, 10).expect("fit should succeed");
    // The clamp `max(0.0)` produces exactly `b = 0`; compare the whole law rather than a bare float.
    assert_eq!(HardApproachLaw { t_hard: 0.0, b: 0.0 }, law);
    Ok(())
  }

  #[test]
  fn test_hard_approach_law_fit_err_without_two_finite_points() -> Result<(), Report> {
    // Only one finite point off the boundary in the fit window: the exponent cannot be regressed.
    let mut y = Array1::from_elem(10, f64::INFINITY);
    y[0] = 0.0;
    let grid = GridFn::from_range_values((0.1, 1.0), y)?;
    let law = HardApproachLaw::fit(&grid, 0.0, Side::Left, 5);
    assert_error!(
      law,
      "Hard-boundary power-law fit on the Left side needs at least two finite grid points off the boundary t_hard=0, found 1"
    );
    Ok(())
  }

  // --- HardApproachLaw::mass (edge-relative) ---
  // Oracle: analytic integral of the power-law density over the gap,
  //   integral_{t_hard}^{t_edge} p_edge * (u/dt_edge)^b du = p_edge * dt_edge / (b + 1).

  #[rustfmt::skip]
  #[rstest]
  // b=1, y_edge=0 (p_edge=1), t_edge=1: integral_0^1 t dt = 1/2
  #[case::power_b1((0.0, 1.0), 0.0,          1.0,                       0.5)]
  // b=2, y_edge=0, t_edge=1: integral_0^1 t^2 dt = 1/3
  #[case::power_b2((0.0, 2.0), 0.0,          1.0,                 1.0 / 3.0)]
  // Flat b=0, y_edge=0, t_edge=0.5: constant p=1 over [0, 0.5] = 0.5
  #[case::flat(    (0.0, 0.0), 0.0,          0.5,                       0.5)]
  // Edge ordinate scales mass: b=1, y_edge=ln2 (p_edge=0.5), t_edge=1: 0.5 * 1 / 2 = 0.25
  #[case::anchored((0.0, 1.0), 2.0_f64.ln(), 1.0,                      0.25)]
  #[trace]
  fn test_hard_approach_law_mass(
    #[case] (t_hard, b): (f64, f64),
    #[case] y_edge: f64,
    #[case] t_edge: f64,
    #[case] expected: f64,
  ) {
    let law = HardApproachLaw { t_hard, b };
    assert_abs_diff_eq!(expected, law.mass(GridEdge { t: t_edge, y: y_edge }), epsilon = 1e-14);
  }

  /// The closed-form mass matches a fine trapezoidal quadrature of the power-law density, an
  /// independent numerical oracle (quadrature vs closed form) across several exponents.
  #[rustfmt::skip]
  #[rstest]
  #[case::b0(  0.0)]
  #[case::b05( 0.5)]
  #[case::b1(  1.0)]
  #[case::b2(  2.0)]
  #[case::b3(  3.0)]
  #[trace]
  fn test_hard_approach_law_mass_matches_quadrature(#[case] b: f64) {
    const N: usize = 2_000_001;
    let (t_hard, t_edge, y_edge) = (0.0, 1.3, 0.4);
    let law = HardApproachLaw { t_hard, b };
    let closed_form = law.mass(GridEdge { t: t_edge, y: y_edge });

    // Trapezoidal integral of the density p(u) = exp(-y(u)) over (t_hard, t_edge], where the density
    // is defined independently of the mass formula from the edge-relative power law.
    let edge = GridEdge { t: t_edge, y: y_edge };
    let dt = (t_edge - t_hard) / (N as f64 - 1.0);
    let quad: f64 = (0..N)
      .map(|i| {
        let u = t_hard + i as f64 * dt;
        let weight = if i == 0 || i == N - 1 { 0.5 } else { 1.0 };
        // Density at u; the boundary point (u == t_hard) has density 0 for b > 0 and p_edge for b = 0.
        let density = (-law.eval(edge, u)).exp();
        weight * density * dt
      })
      .sum();

    assert_abs_diff_eq!(closed_form, quad, epsilon = 1e-6);
  }

  // --- HardApproachLaw::scale ---
  // Scaling every neg-log ordinate by `factor` raises probability to a power, so the exponent scales.

  #[rustfmt::skip]
  #[rstest]
  #[case::double( (0.0, 1.5), 2.0, (0.0, 3.0))]
  #[case::triple( (5.0, 1.0), 3.0, (5.0, 3.0))]
  #[case::flat(   (0.0, 0.0), 4.0, (0.0, 0.0))]
  #[trace]
  fn test_hard_approach_law_scale(
    #[case] (t_hard, b): (f64, f64),
    #[case] factor: f64,
    #[case] (expected_t_hard, expected_b): (f64, f64),
  ) {
    let law = HardApproachLaw { t_hard, b };
    let expected = HardApproachLaw {
      t_hard: expected_t_hard,
      b: expected_b,
    };
    assert_eq!(expected, law.scale(factor));
  }

  // --- HardApproachLaw::negate_arg ---
  // The boundary moves to -t_hard; the exponent is unchanged.

  #[rustfmt::skip]
  #[rstest]
  #[case::divergent((5.0, 1.0),  (-5.0, 1.0))]
  #[case::flat(     (5.0, 0.0),  (-5.0, 0.0))]
  #[case::at_zero(  (0.0, 2.0),  (0.0, 2.0))]
  #[trace]
  fn test_hard_approach_law_negate_arg(#[case] (t_hard, b): (f64, f64), #[case] (expected_t_hard, expected_b): (f64, f64)) {
    let law = HardApproachLaw { t_hard, b };
    let expected = HardApproachLaw {
      t_hard: expected_t_hard,
      b: expected_b,
    };
    assert_eq!(expected, law.negate_arg());
  }

  // --- GridFn integration: approach law in extrapolation ---

  #[test]
  fn test_gridfn_hard_approach_power_law_left() -> Result<(), Report> {
    // Grid [1.0, 5.0], left hard boundary at t=0 with b=1 power-law approach. Grid ordinates follow
    // y(t) = -ln(t): y(1)=0, y(2)=-ln2, etc. The edge-relative law reads the live edge y[0]=0 at
    // t_edge=1, so y(0.5) = 0 - 1*ln(0.5/1) = ln2.
    let grid = GridFn::from_range_values(
      (1.0, 5.0),
      ndarray::array![0.0, -(2.0_f64.ln()), -(3.0_f64.ln()), -(4.0_f64.ln()), -(5.0_f64.ln())],
    )?
    .with_left_extrap(BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 1.0 }));

    // In the gap: y(0.5) = -ln(0.5) = ln(2)
    assert_abs_diff_eq!(2.0_f64.ln(), grid.interp(0.5)?, epsilon = 1e-14);
    // At the hard boundary: +inf (b > 0 means density is zero)
    assert!(grid.interp(0.0)?.is_infinite());
    // Beyond the hard boundary: zero sentinel
    assert_abs_diff_eq!(0.0, grid.interp(-0.5)?, epsilon = 1e-14);
    // On the grid: ordinary interpolation
    assert_abs_diff_eq!(0.0, grid.interp(1.0)?, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gridfn_hard_approach_flat_preserves_boundary() -> Result<(), Report> {
    // b = 0: the law is flat at the live edge ordinate across the whole gap and the boundary.
    let grid = GridFn::from_range_values((0.1, 1.0), ndarray::array![5.0, 5.0, 5.0, 5.0, 5.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 0.0 }));

    assert_abs_diff_eq!(5.0, grid.interp(0.05)?, epsilon = 1e-14);
    assert_abs_diff_eq!(5.0, grid.interp(0.0)?, epsilon = 1e-14);
    assert_abs_diff_eq!(0.0, grid.interp(-0.1)?, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gridfn_shift_y_preserves_gap_shape() -> Result<(), Report> {
    // Edge-relative invariance: a vertical shift of every ordinate shifts the extrapolated gap value
    // by the same delta, with the law object unchanged (it reads the shifted edge).
    let grid = GridFn::from_range_values(
      (1.0, 5.0),
      ndarray::array![0.0, -(2.0_f64.ln()), -(3.0_f64.ln()), -(4.0_f64.ln()), -(5.0_f64.ln())],
    )?
    .with_left_extrap(BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 1.0 }));

    let before = grid.interp(0.5)?;
    let shifted = grid.shift_y(10.0);
    assert_eq!(
      BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 1.0 }),
      shifted.left_extrap()
    );
    assert_abs_diff_eq!(before + 10.0, shifted.interp(0.5)?, epsilon = 1e-13);
    Ok(())
  }

  #[test]
  fn test_gridfn_scale_y_scales_approach_law() -> Result<(), Report> {
    // Scaling every neg-log ordinate by a factor raises probability to a power, so the exponent scales
    // by the same factor while the boundary is unchanged.
    let law = HardApproachLaw { t_hard: 0.0, b: 1.5 };
    let grid = GridFn::from_range_values((1.0, 3.0), ndarray::array![1.0, 2.0, 3.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(law));

    let scaled = grid.scale_y(3.0);
    let scaled_law = scaled
      .left_extrap()
      .approach_law()
      .expect("approach law should be preserved");
    assert_eq!(HardApproachLaw { t_hard: 0.0, b: 4.5 }, scaled_law);
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_clears_approach_law() -> Result<(), Report> {
    let law = HardApproachLaw { t_hard: 0.0, b: 1.0 };
    let grid = GridFn::from_range_values((1.0, 3.0), ndarray::array![1.0, 2.0, 3.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(law));

    let mapped = grid.mapv(|v| v * v);
    assert_eq!(None, mapped.left_extrap().approach_law());
    Ok(())
  }

  #[test]
  fn test_gridfn_negate_arg_swaps_approach_laws() -> Result<(), Report> {
    let left_law = HardApproachLaw { t_hard: 0.0, b: 1.0 };
    let right_law = HardApproachLaw { t_hard: 10.0, b: 2.0 };
    let grid = GridFn::from_range_values((1.0, 9.0), ndarray::array![1.0, 5.0, 9.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(left_law))
      .with_right_extrap(BoundaryBehavior::HardApproach(right_law));

    let negated = grid.negate_arg()?;
    let new_left = negated.left_extrap().approach_law().expect("should have left approach");
    let new_right = negated
      .right_extrap()
      .approach_law()
      .expect("should have right approach");

    // Right law (t_hard=10) reflects to the left: t_hard=-10, exponent unchanged.
    assert_eq!(HardApproachLaw { t_hard: -10.0, b: 2.0 }, new_left);
    // Left law (t_hard=0) reflects to the right: t_hard=0, exponent unchanged.
    assert_eq!(HardApproachLaw { t_hard: 0.0, b: 1.0 }, new_right);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_preserves_approach_law() -> Result<(), Report> {
    let law = HardApproachLaw { t_hard: 0.0, b: 1.0 };
    let grid = GridFn::from_range_values((1.0, 5.0), ndarray::array![2.0, 4.0, 6.0, 8.0, 10.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(law));

    let resampled = grid.resample_range_dx((1.0, 5.0), 0.5)?;
    let resampled_law = resampled
      .left_extrap()
      .approach_law()
      .expect("approach law should be preserved");
    assert_eq!(law, resampled_law);
    Ok(())
  }
}
