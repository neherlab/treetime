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

  // --- HardApproachLaw::eval (edge-relative, per-regime) ---
  // Oracle: `y(t) = y_edge - b * ln(|t - t_hard| / |t_edge - t_hard|) + slope * (t - t_edge)`,
  // evaluated by hand. `from_terms` classifies each `(b, slope)` pair into a regime: `b > 0` with a
  // non-zero slope is `Combined`, `b > 0` alone is `Divergent` (value `+inf` at `t == t_hard`), and
  // `b == 0` is `Finite` (the line `y_edge + slope * (t - t_edge)`, finite at the boundary).

  #[rustfmt::skip]
  #[rstest]
  // b = 1, slope = 0 -> Divergent, boundary at 0, edge (t_edge=1, y_edge=0): y(t) = -ln(t)
  #[case::b1_edge(       (0.0, 1.0, 0.0), (0.0, 1.0), 1.0,                    0.0)]
  #[case::b1_half(       (0.0, 1.0, 0.0), (0.0, 1.0), 0.5,          2.0_f64.ln())]
  #[case::b1_boundary(   (0.0, 1.0, 0.0), (0.0, 1.0), 0.0,          f64::INFINITY)]
  // b = 2, slope = 0 -> Divergent, boundary at 0, edge (1, 0): y(t) = -2 ln(t)
  #[case::b2_half(       (0.0, 2.0, 0.0), (0.0, 1.0), 0.5,    2.0 * 2.0_f64.ln())]
  // Non-zero edge ordinate anchors the law: b=1, edge (t_edge=2, y_edge=3), y(1) = 3 - ln(1/2)
  #[case::anchored(      (0.0, 1.0, 0.0), (3.0, 2.0), 1.0,      3.0 + 2.0_f64.ln())]
  // b = 0, slope = 2 -> Finite line y = y_edge + slope*(t - t_edge), edge (t_edge=1, y_edge=2.5)
  #[case::line_interior( (0.0, 0.0, 2.0), (2.5, 1.0), 0.3,                    1.1)]
  #[case::line_boundary( (0.0, 0.0, 2.0), (2.5, 1.0), 0.0,                    0.5)]
  // b = 0, slope = 0 -> Finite (flat) at the edge ordinate, everywhere including the boundary
  #[case::flat_gap(      (0.0, 0.0, 0.0), (2.5, 1.0), 0.3,                    2.5)]
  #[case::flat_boundary( (0.0, 0.0, 0.0), (2.5, 1.0), 0.0,                    2.5)]
  // b=1, slope=1 -> Combined (product of divergent and finite messages), edge (1, 0):
  // y(0.5) = -ln(0.5) + 1*(0.5 - 1) = ln2 - 0.5
  #[case::composed(      (0.0, 1.0, 1.0), (0.0, 1.0), 0.5,       2.0_f64.ln() - 0.5)]
  // Combined law still diverges at the boundary
  #[case::composed_boundary((0.0, 1.0, 1.0), (0.0, 1.0), 0.0,    f64::INFINITY)]
  // Non-zero t_hard, b = 1, slope = 0 -> Divergent, edge (t_edge=6, y_edge=1): y(t) = 1 - ln(|t-5|/1)
  #[case::offset_edge(   (5.0, 1.0, 0.0), (1.0, 6.0), 6.0,                    1.0)]
  #[case::offset_boundary((5.0, 1.0, 0.0),(1.0, 6.0), 5.0,          f64::INFINITY)]
  #[trace]
  fn test_hard_approach_law_eval(
    #[case] (t_hard, b, slope): (f64, f64, f64),
    #[case] (y_edge, t_edge): (f64, f64),
    #[case] t: f64,
    #[case] expected: f64,
  ) {
    let law = HardApproachLaw::from_terms(t_hard, b, slope);
    let actual = law.eval(y_edge, t_edge, t);
    if expected.is_infinite() {
      assert!(actual.is_infinite() && actual > 0.0, "expected +inf, got {actual}");
    } else {
      assert_abs_diff_eq!(expected, actual, epsilon = 1e-14);
    }
  }

  // --- HardApproachLaw::fit_linear / fit_log_power_law ---
  // `fit_log_power_law` fits an exponent `b` from the grid (divergent boundary); `fit_linear(y_hard)`
  // is the exact line from the boundary to the grid edge (`Finite`, no regression).

  // Divergent case: exact neg-log power-law data, fit must recover the exponent `b` (the discarded
  // intercept `a` does not affect the edge-relative law) and land in the `Divergent` regime.
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
    let law = HardApproachLaw::fit_log_power_law(&grid, t_hard, Side::Left, 10).expect("fit should succeed");
    assert_abs_diff_eq!(t_hard, law.t_hard, epsilon = 1e-14);
    let Approach::Divergent { b: fitted } = law.shape else {
      panic!("expected Divergent regime, got {:?}", law.shape);
    };
    assert_abs_diff_eq!(b, fitted, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_hard_approach_law_fit_right_side() -> Result<(), Report> {
    // y(t) = 2.0 - 1.5 * ln(10 - t), divergent right boundary at t_hard = 10.0 (y_hard = +inf), fit
    // from rightmost points.
    let t_hard = 10.0;
    let b = 1.5;
    let y = Array1::linspace(8.0, 9.9, 20).mapv(|t: f64| 2.0 - b * (t_hard - t).ln());
    let grid = GridFn::from_range_values((8.0, 9.9), y)?;
    let law = HardApproachLaw::fit_log_power_law(&grid, t_hard, Side::Right, 10).expect("fit should succeed");
    let Approach::Divergent { b: fitted } = law.shape else {
      panic!("expected Divergent regime, got {:?}", law.shape);
    };
    assert_abs_diff_eq!(b, fitted, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_hard_approach_law_fit_finite_recovers_slope() -> Result<(), Report> {
    // A zero-mutation branch has a finite, maximal density at the boundary. The neg-log curve is the
    // line `y = 0.5 t + 1`, so the boundary value is `y_hard = 1.0`. The slope is the exact line from
    // the boundary to the grid edge: (1.05 - 1.0) / (0.1 - 0) = 0.5.
    let y = Array1::linspace(0.1, 1.0, 20).mapv(|t: f64| 0.5 * t + 1.0);
    let grid = GridFn::from_range_values((0.1, 1.0), y)?;
    let law = HardApproachLaw::fit_linear(&grid, 0.0, Side::Left, 1.0).expect("fit should succeed");
    let Approach::Finite { slope: fitted } = law.shape else {
      panic!("expected Finite regime, got {:?}", law.shape);
    };
    assert_abs_diff_eq!(0.5, fitted, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_hard_approach_law_fit_finite_anchors_line_at_boundary() -> Result<(), Report> {
    // The finite case draws the exact line between the boundary `(t_hard, y_hard)` and the grid
    // edge, independent of the interior points. Edge is y = 2.0 at t = 0.2, y_hard = 0.5, so the
    // slope is (2.0 - 0.5) / (0.2 - 0) = 7.5. `eval` at the boundary returns y_hard exactly.
    let grid = GridFn::from_range_values((0.2, 1.0), ndarray::array![2.0, 3.0, 4.0, 5.0, 6.0])?;
    let law = HardApproachLaw::fit_linear(&grid, 0.0, Side::Left, 0.5).expect("fit should succeed");
    assert_eq!(Approach::Finite { slope: 7.5 }, law.shape);
    assert_abs_diff_eq!(0.5, law.eval(2.0, 0.2, 0.0), epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_hard_approach_law_fit_err_divergent_without_finite_points() -> Result<(), Report> {
    // Divergent boundary (y_hard = +inf): the left edge is finite but every interior point is +inf,
    // so the regression window has only one finite point and cannot fit an exponent.
    let mut y = Array1::from_elem(10, f64::INFINITY);
    y[0] = 0.0;
    let grid = GridFn::from_range_values((0.1, 1.0), y)?;
    let law = HardApproachLaw::fit_log_power_law(&grid, 0.0, Side::Left, 5);
    assert_error!(
      law,
      "Hard-boundary power-law fit on the Left side needs at least two finite grid points off the boundary t_hard=0, found 1"
    );
    Ok(())
  }

  #[test]
  fn test_hard_approach_law_fit_err_finite_boundary_with_nonfinite_edge() -> Result<(), Report> {
    // Finite boundary (y_hard = 0.5) but a non-finite edge ordinate cannot anchor a line.
    let grid = GridFn::from_range_values((0.1, 1.0), Array1::from_elem(10, f64::INFINITY))?;
    let law = HardApproachLaw::fit_linear(&grid, 0.0, Side::Left, 0.5);
    assert_error!(
      law,
      "Grid Left edge ordinate is non-finite (y=inf at t=0.1); cannot anchor an edge-relative law"
    );
    Ok(())
  }

  // --- HardApproachLaw::mass (edge-relative) ---
  // Oracle: analytic integral of exp(-y(t)) over the gap.
  //  - Divergent (b > 0):            p_edge * |t_edge - t_hard| / (b + 1)
  //  - Finite, slope != 0:           p_edge * (exp(slope * dt_edge) - 1) / slope
  //  - Finite, slope == 0:           p_edge * dt_edge
  //  - Combined (b > 0, slope != 0): p_edge * |t_edge - t_hard| / (b + 1), slope term dropped

  #[rustfmt::skip]
  #[rstest]
  // b=1, slope=0 -> Divergent, y_edge=0 (p_edge=1), t_edge=1: integral_0^1 t dt = 1/2
  #[case::power_b1((0.0, 1.0, 0.0), 0.0,          1.0,                             0.5)]
  // b=2, slope=0 -> Divergent, y_edge=0, t_edge=1: integral_0^1 t^2 dt = 1/3
  #[case::power_b2((0.0, 2.0, 0.0), 0.0,          1.0,                       1.0 / 3.0)]
  // b=0, slope=0 -> Finite (flat), y_edge=0, t_edge=0.5: constant p=1 over [0, 0.5] = 0.5
  #[case::flat(    (0.0, 0.0, 0.0), 0.0,          0.5,                             0.5)]
  // b=0, slope=2 -> Finite line, y_edge=0 (p_edge=1), t_edge=0.5: (exp(2*0.5) - 1) / 2 = (e - 1)/2
  #[case::line(    (0.0, 0.0, 2.0), 0.0,          0.5,             1.0_f64.exp_m1() / 2.0)]
  // Edge ordinate scales mass: b=1, slope=0, y_edge=ln2 (p_edge=0.5), t_edge=1: 0.5 * 1 / 2 = 0.25
  #[case::anchored((0.0, 1.0, 0.0), 2.0_f64.ln(), 1.0,                            0.25)]
  // b=1, slope=2 -> Combined: mass drops the slope term, so it equals the Divergent b=1 mass (0.5).
  // This locks the current approximation; the slope-inclusive integral is a pending decision.
  #[case::combined_slope_dropped((0.0, 1.0, 2.0), 0.0, 1.0,                        0.5)]
  #[trace]
  fn test_hard_approach_law_mass(
    #[case] (t_hard, b, slope): (f64, f64, f64),
    #[case] y_edge: f64,
    #[case] t_edge: f64,
    #[case] expected: f64,
  ) {
    let law = HardApproachLaw::from_terms(t_hard, b, slope);
    assert_abs_diff_eq!(expected, law.mass(y_edge, t_edge), epsilon = 1e-14);
  }

  // --- HardApproachLaw::compose_multiply ---

  #[test]
  fn test_hard_approach_law_compose_multiply_shape_params_add() {
    // A product of a divergent (b>0) and a finite (slope>0) message carries both terms and lands in
    // the Combined regime with the terms summed.
    let a = HardApproachLaw::from_terms(0.0, 1.0, 0.5);
    let b = HardApproachLaw::from_terms(0.0, 2.0, 1.5);
    let expected = HardApproachLaw {
      t_hard: 0.0,
      shape: Approach::Combined { b: 3.0, slope: 2.0 },
    };
    assert_eq!(expected, a.compose_multiply(&b));
  }

  // --- HardApproachLaw::negate_arg ---

  #[test]
  fn test_hard_approach_law_negate_arg() {
    let law = HardApproachLaw::from_terms(5.0, 1.0, 2.0);
    let expected = HardApproachLaw {
      t_hard: -5.0,
      shape: Approach::Combined { b: 1.0, slope: -2.0 },
    };
    assert_eq!(expected, law.negate_arg());
  }

  // --- GridFn integration: approach law in extrapolation ---

  #[test]
  fn test_gridfn_hard_approach_power_law_left() -> Result<(), Report> {
    // Grid [1.0, 5.0], left hard boundary at t=0 with b=1 power-law (Divergent) approach.
    // Grid ordinates follow y(t) = -ln(t): y(1)=0, y(2)=-ln2, etc. The edge-relative law reads the
    // live edge y[0]=0 at t_edge=1, so y(0.5) = 0 - 1*ln(0.5/1) = ln2.
    let grid = GridFn::from_range_values(
      (1.0, 5.0),
      ndarray::array![0.0, -(2.0_f64.ln()), -(3.0_f64.ln()), -(4.0_f64.ln()), -(5.0_f64.ln())],
    )?
    .with_left_extrap(BoundaryBehavior::HardApproach(HardApproachLaw {
      t_hard: 0.0,
      shape: Approach::Divergent { b: 1.0 },
    }));

    // In the gap: y(0.5) = -ln(0.5) = ln(2)
    assert_abs_diff_eq!(2.0_f64.ln(), grid.interp(0.5)?, epsilon = 1e-14);
    // At the hard boundary: +inf (Divergent means density is zero)
    assert!(grid.interp(0.0)?.is_infinite());
    // Beyond the hard boundary: zero sentinel
    assert_abs_diff_eq!(0.0, grid.interp(-0.5)?, epsilon = 1e-14);
    // On the grid: ordinary interpolation
    assert_abs_diff_eq!(0.0, grid.interp(1.0)?, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gridfn_hard_approach_finite_line_reaches_boundary() -> Result<(), Report> {
    // Finite, slope > 0: the line carries the mode to the boundary. Grid [1.0, 3.0] with a
    // live edge y[0] = 4.0 at t_edge = 1.0 and slope = 2.0. In the gap the law is the line
    // y = 4.0 + 2.0*(t - 1.0); at the boundary t = 0 it reaches the finite value 2.0 (not y_edge).
    let grid = GridFn::from_range_values((1.0, 3.0), ndarray::array![4.0, 6.0, 8.0])?.with_left_extrap(
      BoundaryBehavior::HardApproach(HardApproachLaw {
        t_hard: 0.0,
        shape: Approach::Finite { slope: 2.0 },
      }),
    );

    assert_abs_diff_eq!(3.0, grid.interp(0.5)?, epsilon = 1e-14);
    assert_abs_diff_eq!(2.0, grid.interp(0.0)?, epsilon = 1e-14);
    assert_abs_diff_eq!(0.0, grid.interp(-0.1)?, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gridfn_hard_approach_flat_preserves_boundary() -> Result<(), Report> {
    // Finite, slope = 0: the law is flat at the live edge ordinate across the whole gap and boundary.
    let grid = GridFn::from_range_values((0.1, 1.0), ndarray::array![5.0, 5.0, 5.0, 5.0, 5.0])?.with_left_extrap(
      BoundaryBehavior::HardApproach(HardApproachLaw {
        t_hard: 0.0,
        shape: Approach::Finite { slope: 0.0 },
      }),
    );

    assert_abs_diff_eq!(5.0, grid.interp(0.05)?, epsilon = 1e-14);
    assert_abs_diff_eq!(5.0, grid.interp(0.0)?, epsilon = 1e-14);
    assert_abs_diff_eq!(0.0, grid.interp(-0.1)?, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gridfn_shift_y_preserves_gap_shape() -> Result<(), Report> {
    // Edge-relative invariance: a vertical shift of every ordinate shifts the extrapolated gap
    // value by the same delta, with the law object unchanged (it reads the shifted edge).
    let grid = GridFn::from_range_values(
      (1.0, 5.0),
      ndarray::array![0.0, -(2.0_f64.ln()), -(3.0_f64.ln()), -(4.0_f64.ln()), -(5.0_f64.ln())],
    )?
    .with_left_extrap(BoundaryBehavior::HardApproach(HardApproachLaw {
      t_hard: 0.0,
      shape: Approach::Divergent { b: 1.0 },
    }));

    let before = grid.interp(0.5)?;
    let shifted = grid.shift_y(10.0);
    assert_eq!(
      BoundaryBehavior::HardApproach(HardApproachLaw {
        t_hard: 0.0,
        shape: Approach::Divergent { b: 1.0 },
      }),
      shifted.left_extrap()
    );
    assert_abs_diff_eq!(before + 10.0, shifted.interp(0.5)?, epsilon = 1e-13);
    Ok(())
  }

  #[test]
  fn test_gridfn_scale_y_scales_approach_law() -> Result<(), Report> {
    let law = HardApproachLaw {
      t_hard: 0.0,
      shape: Approach::Combined { b: 1.5, slope: 2.0 },
    };
    let grid = GridFn::from_range_values((1.0, 3.0), ndarray::array![1.0, 2.0, 3.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(law));

    let scaled = grid.scale_y(3.0);
    let scaled_law = scaled
      .left_extrap()
      .approach_law()
      .expect("approach law should be preserved");
    // Scaling ordinates raises probability to a power, so both shape terms scale by the factor.
    assert_eq!(
      HardApproachLaw {
        t_hard: 0.0,
        shape: Approach::Combined { b: 4.5, slope: 6.0 },
      },
      scaled_law
    );
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_clears_approach_law() -> Result<(), Report> {
    let law = HardApproachLaw {
      t_hard: 0.0,
      shape: Approach::Divergent { b: 1.0 },
    };
    let grid = GridFn::from_range_values((1.0, 3.0), ndarray::array![1.0, 2.0, 3.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(law));

    let mapped = grid.mapv(|v| v * v);
    assert_eq!(None, mapped.left_extrap().approach_law());
    Ok(())
  }

  #[test]
  fn test_gridfn_negate_arg_swaps_approach_laws() -> Result<(), Report> {
    let left_law = HardApproachLaw {
      t_hard: 0.0,
      shape: Approach::Combined { b: 1.0, slope: 2.0 },
    };
    let right_law = HardApproachLaw {
      t_hard: 10.0,
      shape: Approach::Combined { b: 2.0, slope: -1.0 },
    };
    let grid = GridFn::from_range_values((1.0, 9.0), ndarray::array![1.0, 5.0, 9.0])?
      .with_left_extrap(BoundaryBehavior::HardApproach(left_law))
      .with_right_extrap(BoundaryBehavior::HardApproach(right_law));

    let negated = grid.negate_arg()?;
    let new_left = negated.left_extrap().approach_law().expect("should have left approach");
    let new_right = negated
      .right_extrap()
      .approach_law()
      .expect("should have right approach");

    // Right law (t_hard=10) reflects to the left: t_hard=-10, exponent unchanged, slope flips.
    assert_eq!(
      HardApproachLaw {
        t_hard: -10.0,
        shape: Approach::Combined { b: 2.0, slope: 1.0 },
      },
      new_left
    );
    // Left law (t_hard=0) reflects to the right: t_hard=0, exponent unchanged, slope flips.
    assert_eq!(
      HardApproachLaw {
        t_hard: 0.0,
        shape: Approach::Combined { b: 1.0, slope: -2.0 },
      },
      new_right
    );
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_preserves_approach_law() -> Result<(), Report> {
    let law = HardApproachLaw {
      t_hard: 0.0,
      shape: Approach::Combined { b: 1.0, slope: 0.5 },
    };
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
