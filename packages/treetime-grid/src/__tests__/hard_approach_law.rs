#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::*;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  fn make_neglog_power_law_grid(t_hard: f64, a: f64, b: f64, x_min: f64, x_max: f64, n: usize) -> GridFn<f64> {
    let y = Array1::linspace(x_min, x_max, n).mapv(|t| a - b * (t - t_hard).abs().ln());
    GridFn::from_range_values((x_min, x_max), y).unwrap()
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::b1_edge(       (0.0, 1.0), (0.0, 1.0), 1.0,                    0.0)]
  #[case::b1_half(       (0.0, 1.0), (0.0, 1.0), 0.5,          2.0_f64.ln())]
  #[case::b1_boundary(   (0.0, 1.0), (0.0, 1.0), 0.0,          f64::INFINITY)]
  #[case::b2_half(       (0.0, 2.0), (0.0, 1.0), 0.5,    2.0 * 2.0_f64.ln())]
  #[case::anchored(      (0.0, 1.0), (3.0, 2.0), 1.0,      3.0 + 2.0_f64.ln())]
  #[case::flat_gap(      (0.0, 0.0), (2.5, 1.0), 0.3,                    2.5)]
  #[case::flat_boundary( (0.0, 0.0), (2.5, 1.0), 0.0,                    2.5)]
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
    let y = Array1::linspace(0.1, 1.0, 20).mapv(|t: f64| 2.0 + (t - 0.0).abs().ln());
    let grid = GridFn::from_range_values((0.1, 1.0), y)?;
    let law = HardApproachLaw::fit(&grid, 0.0, Side::Left, 10).expect("fit should succeed");
    assert_eq!(HardApproachLaw { t_hard: 0.0, b: 0.0 }, law);
    Ok(())
  }

  #[test]
  fn test_hard_approach_law_fit_err_without_two_finite_points() -> Result<(), Report> {
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

  #[rustfmt::skip]
  #[rstest]
  #[case::power_b1((0.0, 1.0), 0.0,          1.0,                       0.5)]
  #[case::power_b2((0.0, 2.0), 0.0,          1.0,                 1.0 / 3.0)]
  #[case::flat(    (0.0, 0.0), 0.0,          0.5,                       0.5)]
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

    let edge = GridEdge { t: t_edge, y: y_edge };
    let dt = (t_edge - t_hard) / (N as f64 - 1.0);
    let quad: f64 = (0..N)
      .map(|i| {
        let u = t_hard + i as f64 * dt;
        let weight = if i == 0 || i == N - 1 { 0.5 } else { 1.0 };
        let density = (-law.eval(edge, u)).exp();
        weight * density * dt
      })
      .sum();

    assert_abs_diff_eq!(closed_form, quad, epsilon = 1e-6);
  }

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

  #[test]
  fn test_gridfn_hard_approach_power_law_left() -> Result<(), Report> {
    let grid = GridFn::from_range_values(
      (1.0, 5.0),
      ndarray::array![0.0, -(2.0_f64.ln()), -(3.0_f64.ln()), -(4.0_f64.ln()), -(5.0_f64.ln())],
    )?;
    let behavior = BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 1.0 });

    assert_abs_diff_eq!(
      2.0_f64.ln(),
      grid.interp_with_extrap(0.5, behavior, BoundaryBehavior::Error)?,
      epsilon = 1e-14
    );
    assert!(
      grid
        .interp_with_extrap(0.0, behavior, BoundaryBehavior::Error)?
        .is_infinite()
    );
    assert_abs_diff_eq!(
      0.0,
      grid.interp_with_extrap(-0.5, behavior, BoundaryBehavior::Error)?,
      epsilon = 1e-14
    );
    assert_abs_diff_eq!(0.0, grid.interp(1.0)?, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gridfn_hard_approach_flat_preserves_boundary() -> Result<(), Report> {
    let grid = GridFn::from_range_values((0.1, 1.0), ndarray::array![5.0, 5.0, 5.0, 5.0, 5.0])?;
    let behavior = BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 0.0 });

    assert_abs_diff_eq!(
      5.0,
      grid.interp_with_extrap(0.05, behavior, BoundaryBehavior::Error)?,
      epsilon = 1e-14
    );
    assert_abs_diff_eq!(
      5.0,
      grid.interp_with_extrap(0.0, behavior, BoundaryBehavior::Error)?,
      epsilon = 1e-14
    );
    assert_abs_diff_eq!(
      0.0,
      grid.interp_with_extrap(-0.1, behavior, BoundaryBehavior::Error)?,
      epsilon = 1e-14
    );
    Ok(())
  }

  #[test]
  fn test_gridfn_shift_y_preserves_gap_shape() -> Result<(), Report> {
    let grid = GridFn::from_range_values(
      (1.0, 5.0),
      ndarray::array![0.0, -(2.0_f64.ln()), -(3.0_f64.ln()), -(4.0_f64.ln()), -(5.0_f64.ln())],
    )?;
    let behavior = BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 1.0 });

    let before = grid.interp_with_extrap(0.5, behavior, BoundaryBehavior::Error)?;
    let shifted = grid.shift_y(10.0);
    assert_abs_diff_eq!(
      before + 10.0,
      shifted.interp_with_extrap(0.5, behavior, BoundaryBehavior::Error)?,
      epsilon = 1e-13
    );
    Ok(())
  }
}
