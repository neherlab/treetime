#[cfg(test)]
mod tests {
  use crate::*;
  use approx::{assert_abs_diff_eq, assert_ulps_eq};
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  /// Build a GridFn whose stored ordinates follow neg-log of a power law:
  /// `y(t) = a - b * ln|t - t_hard|`, for testing that the fit recovers the known parameters.
  fn make_neglog_power_law_grid(t_hard: f64, a: f64, b: f64, x_min: f64, x_max: f64, n: usize) -> GridFn<f64> {
    let y = Array1::linspace(x_min, x_max, n).mapv(|t| a - b * (t - t_hard).abs().ln());
    GridFn::from_range_values((x_min, x_max), y).unwrap()
  }

  /// Build a GridFn whose stored ordinates follow a straight line:
  /// `y(t) = slope * t + intercept`, for testing the linear (n=0) case.
  fn make_linear_grid(slope: f64, intercept: f64, x_min: f64, x_max: f64, n: usize) -> GridFn<f64> {
    let y = Array1::linspace(x_min, x_max, n).mapv(|t| slope * t + intercept);
    GridFn::from_range_values((x_min, x_max), y).unwrap()
  }

  // --- HardApproachLaw::eval ---
  // Oracle: `y(t) = a - b * ln|t - t_hard| + slope * (t - t_hard)`, evaluated by hand.

  #[rustfmt::skip]
  #[rstest]
  // b=0, slope=0: constant (degenerate flat case)
  #[case::flat_at_boundary(  (0.0, 2.5, 0.0, 0.0), 0.0,   2.5)]
  #[case::flat_away(         (0.0, 2.5, 0.0, 0.0), 0.5,   2.5)]
  // b=0, slope>0: linear (zero-mutation branch, neg-log density y = a + slope*t)
  #[case::linear_at_hard(    (0.0, 1.0, 0.0, 3.0), 0.0,   1.0)]
  #[case::linear_unit(       (0.0, 1.0, 0.0, 3.0), 1.0,   4.0)]
  #[case::linear_half(       (0.0, 1.0, 0.0, 3.0), 0.5,   2.5)]
  // b=1, slope=0: one-mutation power law, neg-log y = a - ln|dt|
  #[case::b1_at_boundary(    (0.0, 0.0, 1.0, 0.0), 0.0,   f64::INFINITY)]
  #[case::b1_unit(           (0.0, 0.0, 1.0, 0.0), 1.0,   0.0)]
  #[case::b1_e(              (0.0, 0.0, 1.0, 0.0), std::f64::consts::E, -1.0)]
  // b=2, slope=0: two-mutation power law, neg-log y = a - 2*ln|dt|
  #[case::b2_at_boundary(    (0.0, 0.0, 2.0, 0.0), 0.0,   f64::INFINITY)]
  #[case::b2_unit(           (0.0, 0.0, 2.0, 0.0), 1.0,   0.0)]
  #[case::b2_half(           (0.0, 0.0, 2.0, 0.0), 0.5,   2.0 * 2.0_f64.ln())]
  // Non-zero t_hard with b > 0
  #[case::offset_at_hard(    (5.0, 1.0, 1.0, 0.0), 5.0,   f64::INFINITY)]
  #[case::offset_away(       (5.0, 1.0, 1.0, 0.0), 6.0,   1.0)]
  // Combined b > 0 and slope > 0 (from composition of zero- and nonzero-mutation messages)
  #[case::combined_unit(     (0.0, 1.0, 1.0, 2.0), 1.0,   1.0 + 2.0)]
  #[trace]
  fn test_hard_approach_law_eval(
    #[case] (t_hard, a, b, slope): (f64, f64, f64, f64),
    #[case] t: f64,
    #[case] expected: f64,
  ) {
    let law = HardApproachLaw { t_hard, a, b, slope };
    let actual = law.eval(t);
    if expected.is_infinite() {
      assert!(actual.is_infinite() && actual > 0.0, "expected +inf, got {actual}");
    } else {
      assert_abs_diff_eq!(expected, actual, epsilon = 1e-14);
    }
  }

  // --- HardApproachLaw::fit ---

  // Power-law case: exact neg-log power-law data, fit must recover (a, b) with slope = 0.
  // Oracle: stored ordinates follow `y = a - b * ln|dt|` exactly.
  #[rustfmt::skip]
  #[rstest]
  #[case::b1(  0.0, 1.0, 1.0)]
  #[case::b2(  0.0, 0.5, 2.0)]
  #[case::b3(  0.0, 2.0, 3.0)]
  #[case::b05( 0.0, 1.0, 0.5)]
  #[trace]
  fn test_hard_approach_law_fit_recovers_power_law(
    #[case] t_hard: f64,
    #[case] a: f64,
    #[case] b: f64,
  ) -> Result<(), Report> {
    let grid = make_neglog_power_law_grid(t_hard, a, b, 0.1, 1.0, 20);
    let law = HardApproachLaw::fit(&grid, t_hard, Side::Left, 10).expect("fit should succeed");
    assert_abs_diff_eq!(b, law.b, epsilon = 1e-10);
    assert_abs_diff_eq!(a, law.a, epsilon = 1e-10);
    assert_abs_diff_eq!(0.0, law.slope, epsilon = 1e-14);
    assert_ulps_eq!(t_hard, law.t_hard);
    Ok(())
  }

  // Linear case: exact linear data, fit must recover (a = intercept, slope) with b = 0.
  // Oracle: stored ordinates follow `y = slope * t + intercept` exactly.
  #[rustfmt::skip]
  #[rstest]
  #[case::flat(          0.0,  5.0)]
  #[case::rising(        3.0,  0.5)]
  #[case::falling(      -2.0,  4.0)]
  #[case::steep(         7.5, -1.0)]
  #[trace]
  fn test_hard_approach_law_fit_recovers_linear(
    #[case] slope: f64,
    #[case] intercept: f64,
  ) -> Result<(), Report> {
    let grid = make_linear_grid(slope, intercept, 0.1, 1.0, 20);
    let law = HardApproachLaw::fit(&grid, 0.0, Side::Left, 10).expect("fit should succeed");
    assert_abs_diff_eq!(0.0, law.b, epsilon = 1e-14);
    assert_abs_diff_eq!(intercept, law.a, epsilon = 1e-10); // a = slope*0 + intercept
    assert_abs_diff_eq!(slope, law.slope, epsilon = 1e-10);
    assert_ulps_eq!(0.0, law.t_hard);
    Ok(())
  }

  #[test]
  fn test_hard_approach_law_fit_right_side() -> Result<(), Report> {
    // y(t) = 2.0 - 1.5 * ln(10 - t), right boundary at t_hard = 10.0, fit from rightmost points.
    let t_hard = 10.0;
    let a = 2.0;
    let b = 1.5;
    let y = Array1::linspace(8.0, 9.9, 20).mapv(|t: f64| a - b * (t_hard - t).ln());
    let grid = GridFn::from_range_values((8.0, 9.9), y)?;
    let law = HardApproachLaw::fit(&grid, t_hard, Side::Right, 10).expect("fit should succeed");
    assert_abs_diff_eq!(b, law.b, epsilon = 1e-10);
    assert_abs_diff_eq!(a, law.a, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_hard_approach_law_fit_returns_none_for_non_finite() -> Result<(), Report> {
    let grid = GridFn::from_range_values((0.1, 1.0), Array1::from_elem(10, f64::INFINITY))?;
    let law = HardApproachLaw::fit(&grid, 0.0, Side::Left, 5);
    assert_eq!(None, law);
    Ok(())
  }

  // --- HardApproachLaw::mass ---
  // Oracle: analytic integral of exp(-y(t)) over the gap.

  #[test]
  fn test_hard_approach_law_mass_power_law() {
    // b = 1, a = 0, slope = 0: p(t) = t^1, mass = integral_0^1 t dt = 1/2.
    let law = HardApproachLaw {
      t_hard: 0.0,
      a: 0.0,
      b: 1.0,
      slope: 0.0,
    };
    assert_abs_diff_eq!(0.5, law.mass(1.0), epsilon = 1e-14);
  }

  #[test]
  fn test_hard_approach_law_mass_power_law_b2() {
    // b = 2, a = 0, slope = 0: p(t) = t^2, mass = integral_0^1 t^2 dt = 1/3.
    let law = HardApproachLaw {
      t_hard: 0.0,
      a: 0.0,
      b: 2.0,
      slope: 0.0,
    };
    assert_abs_diff_eq!(1.0 / 3.0, law.mass(1.0), epsilon = 1e-14);
  }

  #[test]
  fn test_hard_approach_law_mass_linear() {
    // b = 0, slope = 1, a = 0: integral_0^1 exp(-t) dt = 1 - exp(-1).
    let law = HardApproachLaw {
      t_hard: 0.0,
      a: 0.0,
      b: 0.0,
      slope: 1.0,
    };
    assert_abs_diff_eq!(1.0 - (-1.0_f64).exp(), law.mass(1.0), epsilon = 1e-14);
  }

  #[test]
  fn test_hard_approach_law_mass_flat() {
    // b = 0, slope = 0, a = 0: constant p = 1, mass = dt = 0.5.
    let law = HardApproachLaw {
      t_hard: 0.0,
      a: 0.0,
      b: 0.0,
      slope: 0.0,
    };
    assert_abs_diff_eq!(0.5, law.mass(0.5), epsilon = 1e-14);
  }

  // --- HardApproachLaw::compose_multiply ---

  #[test]
  fn test_hard_approach_law_compose_multiply_all_add() {
    let a = HardApproachLaw {
      t_hard: 0.0,
      a: 1.0,
      b: 1.0,
      slope: 2.0,
    };
    let b = HardApproachLaw {
      t_hard: 0.0,
      a: 3.0,
      b: 2.0,
      slope: 0.5,
    };
    let expected = HardApproachLaw {
      t_hard: 0.0,
      a: 4.0,
      b: 3.0,
      slope: 2.5,
    };
    assert_eq!(expected, a.compose_multiply(&b));
  }

  #[test]
  fn test_hard_approach_law_compose_multiply_linear_times_linear() {
    let a = HardApproachLaw {
      t_hard: 0.0,
      a: 2.0,
      b: 0.0,
      slope: 1.0,
    };
    let b = HardApproachLaw {
      t_hard: 0.0,
      a: 5.0,
      b: 0.0,
      slope: 3.0,
    };
    let result = a.compose_multiply(&b);
    assert_abs_diff_eq!(0.0, result.b, epsilon = 1e-14);
    assert_abs_diff_eq!(7.0, result.a, epsilon = 1e-14);
    assert_abs_diff_eq!(4.0, result.slope, epsilon = 1e-14);
  }

  // --- HardApproachLaw::negate_arg ---

  #[test]
  fn test_hard_approach_law_negate_arg() {
    let law = HardApproachLaw {
      t_hard: 5.0,
      a: 2.0,
      b: 1.0,
      slope: 3.0,
    };
    let expected = HardApproachLaw {
      t_hard: -5.0,
      a: 2.0,
      b: 1.0,
      slope: -3.0,
    };
    assert_eq!(expected, law.negate_arg());
  }

  // --- HardApproachLaw::shift_anchor ---

  #[test]
  fn test_hard_approach_law_shift_anchor() {
    let law = HardApproachLaw {
      t_hard: 0.0,
      a: 5.0,
      b: 1.0,
      slope: 2.0,
    };
    let shifted = law.shift_anchor(-3.0);
    assert_abs_diff_eq!(2.0, shifted.a, epsilon = 1e-14);
    assert_abs_diff_eq!(1.0, shifted.b, epsilon = 1e-14);
    assert_abs_diff_eq!(2.0, shifted.slope, epsilon = 1e-14);
  }

  // --- GridFn integration: approach law in extrapolation ---

  #[test]
  fn test_gridfn_hard_approach_power_law_left() -> Result<(), Report> {
    // Grid [1.0, 5.0], left hard boundary at t=0 with b=1 power-law approach.
    // y(t) = a - b*ln(t) = 0 - 1*ln(t) = -ln(t).
    // Grid ordinates follow the same law: y(1)=0, y(2)=-ln2, etc.
    let grid = GridFn::from_range_values(
      (1.0, 5.0),
      ndarray::array![0.0, -(2.0_f64.ln()), -(3.0_f64.ln()), -(4.0_f64.ln()), -(5.0_f64.ln())],
    )?
    .with_left_extrap(BoundaryBehavior::Hard(Some(HardApproachLaw {
      t_hard: 0.0,
      a: 0.0,
      b: 1.0,
      slope: 0.0,
    })));

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
  fn test_gridfn_hard_approach_linear_left() -> Result<(), Report> {
    // Grid [0.1, 0.5], left hard boundary at t=0, linear approach y = 2*t.
    let grid = GridFn::from_range_values(
      (0.1, 0.5),
      ndarray::array![0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    )?
    .with_left_extrap(BoundaryBehavior::Hard(Some(HardApproachLaw {
      t_hard: 0.0,
      a: 0.0,
      b: 0.0,
      slope: 2.0,
    })));

    // In the gap: y(0.05) = 0 + 2*0.05 = 0.1
    assert_abs_diff_eq!(0.1, grid.interp(0.05)?, epsilon = 1e-14);
    // At the hard boundary: y(0) = 0 (finite, mode preserved)
    assert_abs_diff_eq!(0.0, grid.interp(0.0)?, epsilon = 1e-14);
    // Beyond: zero sentinel
    assert_abs_diff_eq!(0.0, grid.interp(-0.1)?, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gridfn_hard_approach_flat_preserves_boundary() -> Result<(), Report> {
    // b = 0, slope = 0: constant ordinate across the gap.
    let grid = GridFn::from_range_values((0.1, 1.0), ndarray::array![5.0, 5.0, 5.0, 5.0, 5.0])?.with_left_extrap(
      BoundaryBehavior::Hard(Some(HardApproachLaw {
        t_hard: 0.0,
        a: 5.0,
        b: 0.0,
        slope: 0.0,
      })),
    );

    assert_abs_diff_eq!(5.0, grid.interp(0.05)?, epsilon = 1e-14);
    assert_abs_diff_eq!(5.0, grid.interp(0.0)?, epsilon = 1e-14);
    assert_abs_diff_eq!(0.0, grid.interp(-0.1)?, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gridfn_scale_y_scales_approach_law() -> Result<(), Report> {
    let law = HardApproachLaw {
      t_hard: 0.0,
      a: 1.0,
      b: 1.5,
      slope: 2.0,
    };
    let grid = GridFn::from_range_values((1.0, 3.0), ndarray::array![1.0, 2.0, 3.0])?
      .with_left_extrap(BoundaryBehavior::Hard(Some(law)));

    let scaled = grid.scale_y(3.0);
    let scaled_law = scaled
      .left_extrap()
      .approach_law()
      .expect("approach law should be preserved");
    let expected = HardApproachLaw {
      t_hard: 0.0,
      a: 3.0,
      b: 1.5,
      slope: 6.0,
    };
    assert_eq!(expected, scaled_law);
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_clears_approach_law() -> Result<(), Report> {
    let law = HardApproachLaw {
      t_hard: 0.0,
      a: 1.0,
      b: 1.0,
      slope: 0.0,
    };
    let grid = GridFn::from_range_values((1.0, 3.0), ndarray::array![1.0, 2.0, 3.0])?
      .with_left_extrap(BoundaryBehavior::Hard(Some(law)));

    let mapped = grid.mapv(|v| v * v);
    assert_eq!(None, mapped.left_extrap().approach_law());
    Ok(())
  }

  #[test]
  fn test_gridfn_negate_arg_swaps_approach_laws() -> Result<(), Report> {
    let left_law = HardApproachLaw {
      t_hard: 0.0,
      a: 1.0,
      b: 1.0,
      slope: 2.0,
    };
    let right_law = HardApproachLaw {
      t_hard: 10.0,
      a: 3.0,
      b: 2.0,
      slope: -1.0,
    };
    let grid = GridFn::from_range_values((1.0, 9.0), ndarray::array![1.0, 5.0, 9.0])?
      .with_left_extrap(BoundaryBehavior::Hard(Some(left_law)))
      .with_right_extrap(BoundaryBehavior::Hard(Some(right_law)));

    let negated = grid.negate_arg()?;
    let new_left = negated.left_extrap().approach_law().expect("should have left approach");
    let new_right = negated
      .right_extrap()
      .approach_law()
      .expect("should have right approach");

    // Right law (t_hard=10) reflects to the left: t_hard=-10, a unchanged, slope flips.
    assert_eq!(
      HardApproachLaw {
        t_hard: -10.0,
        a: 3.0,
        b: 2.0,
        slope: 1.0
      },
      new_left
    );
    // Left law (t_hard=0) reflects to the right: t_hard=0, a unchanged, slope flips.
    assert_eq!(
      HardApproachLaw {
        t_hard: 0.0,
        a: 1.0,
        b: 1.0,
        slope: -2.0
      },
      new_right
    );
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_preserves_approach_law() -> Result<(), Report> {
    let law = HardApproachLaw {
      t_hard: 0.0,
      a: 1.0,
      b: 1.0,
      slope: 2.0,
    };
    let grid = GridFn::from_range_values((1.0, 5.0), ndarray::array![2.0, 4.0, 6.0, 8.0, 10.0])?
      .with_left_extrap(BoundaryBehavior::Hard(Some(law)));

    let resampled = grid.resample_range_dx((1.0, 5.0), 0.5)?;
    let resampled_law = resampled
      .left_extrap()
      .approach_law()
      .expect("approach law should be preserved");
    assert_eq!(law, resampled_law);
    Ok(())
  }
}
