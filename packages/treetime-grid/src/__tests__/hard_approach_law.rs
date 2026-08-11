#[cfg(test)]
mod tests {
  use crate::*;
  use approx::{assert_abs_diff_eq, assert_ulps_eq};
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  /// Build a GridFn whose y-values follow p(t) = C * (t - t_hard)^b exactly,
  /// for testing that HardApproachLaw::fit recovers the known parameters.
  fn make_power_law_grid(t_hard: f64, coeff: f64, exponent: f64, x_min: f64, x_max: f64, n: usize) -> GridFn<f64> {
    let y = Array1::linspace(x_min, x_max, n).mapv(|t| coeff * (t - t_hard).abs().powf(exponent));
    GridFn::from_range_values((x_min, x_max), y).unwrap()
  }

  // --- HardApproachLaw::eval ---

  #[rustfmt::skip]
  #[rstest]
  // b=0: constant at boundary (zero-mutation case)
  #[case::b0_at_boundary(   0.0, 2.5, 0.0, 0.0,   2.5)]
  #[case::b0_away(          0.0, 2.5, 0.0, 0.5,   2.5)]
  #[case::b0_far(           0.0, 2.5, 0.0, 10.0,  2.5)]
  // b=1: linear decay (one-mutation case)
  #[case::b1_at_boundary(   0.0, 3.0, 1.0, 0.0,   0.0)]
  #[case::b1_unit(          0.0, 3.0, 1.0, 1.0,   3.0)]
  #[case::b1_half(          0.0, 3.0, 1.0, 0.5,   1.5)]
  // b=2: quadratic decay (two-mutation case)
  #[case::b2_at_boundary(   0.0, 1.0, 2.0, 0.0,   0.0)]
  #[case::b2_unit(          0.0, 1.0, 2.0, 1.0,   1.0)]
  #[case::b2_half(          0.0, 1.0, 2.0, 0.5,   0.25)]
  // Non-zero t_hard
  #[case::offset_at_hard(   5.0, 2.0, 1.0, 5.0,   0.0)]
  #[case::offset_away(      5.0, 2.0, 1.0, 6.0,   2.0)]
  #[trace]
  fn test_approach_law_eval(
    #[case] t_hard: f64,
    #[case] coeff: f64,
    #[case] exponent: f64,
    #[case] t: f64,
    #[case] expected: f64,
  ) {
    let law = HardApproachLaw { t_hard, coeff, exponent };
    let actual = law.eval(t);
    assert_abs_diff_eq!(expected, actual, epsilon = 1e-14);
  }

  // --- HardApproachLaw::fit ---

  #[rustfmt::skip]
  #[rstest]
  // Exact power-law data: fit should recover parameters precisely
  #[case::b0_const(    0.0, 5.0, 0.0)]
  #[case::b1_linear(   0.0, 3.0, 1.0)]
  #[case::b2_quadratic(0.0, 1.0, 2.0)]
  #[case::b3_cubic(    0.0, 2.0, 3.0)]
  #[case::b05_sqrt(    0.0, 4.0, 0.5)]
  #[trace]
  fn test_approach_law_fit_recovers_exact_power_law(
    #[case] t_hard: f64,
    #[case] coeff: f64,
    #[case] exponent: f64,
  ) -> Result<(), Report> {
    let grid = make_power_law_grid(t_hard, coeff, exponent, 0.1, 1.0, 20);
    let law = HardApproachLaw::fit(&grid, t_hard, Side::Left, 10).expect("fit should succeed");
    assert_abs_diff_eq!(exponent, law.exponent, epsilon = 1e-10);
    assert_abs_diff_eq!(coeff, law.coeff, epsilon = 1e-10);
    assert_ulps_eq!(t_hard, law.t_hard);
    Ok(())
  }

  #[test]
  fn test_approach_law_fit_right_side() -> Result<(), Report> {
    // p(t) = 2.0 * (t_hard - t)^1.5, t_hard=10.0, grid [8.0, 9.9]
    let t_hard = 10.0;
    let coeff = 2.0;
    let exponent = 1.5;
    let y = Array1::linspace(8.0, 9.9, 20).mapv(|t: f64| coeff * (t_hard - t).powf(exponent));
    let grid = GridFn::from_range_values((8.0, 9.9), y)?;
    let law = HardApproachLaw::fit(&grid, t_hard, Side::Right, 10).expect("fit should succeed");
    assert_abs_diff_eq!(exponent, law.exponent, epsilon = 1e-10);
    assert_abs_diff_eq!(coeff, law.coeff, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_approach_law_fit_clamps_negative_exponent_to_zero() -> Result<(), Report> {
    // Density that increases toward the boundary (like exp(-t) near t=0):
    // log(p) vs log(t) has negative slope, so fit would give b < 0.
    // The approach law clamps to b=0.
    let y = Array1::linspace(0.1, 1.0, 10).mapv(|t: f64| (-t).exp());
    let grid = GridFn::from_range_values((0.1, 1.0), y)?;
    let law = HardApproachLaw::fit(&grid, 0.0, Side::Left, 5).expect("fit should succeed");
    assert_abs_diff_eq!(0.0, law.exponent, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_approach_law_fit_returns_none_for_all_zeros() -> Result<(), Report> {
    let grid = GridFn::from_range_values((0.1, 1.0), Array1::zeros(10))?;
    let law = HardApproachLaw::fit(&grid, 0.0, Side::Left, 5);
    assert_eq!(None, law);
    Ok(())
  }

  // --- HardApproachLaw::mass ---

  #[test]
  fn test_approach_law_mass_b0() {
    // p(t) = C, mass = C * dt
    let law = HardApproachLaw {
      t_hard: 0.0,
      coeff: 3.0,
      exponent: 0.0,
    };
    let mass = law.mass(0.5);
    // 3.0 * 0.5^1 / 1.0 = 1.5
    assert_abs_diff_eq!(1.5, mass, epsilon = 1e-14);
  }

  #[test]
  fn test_approach_law_mass_b1() {
    // p(t) = C * t, mass = C * t^2 / 2
    let law = HardApproachLaw {
      t_hard: 0.0,
      coeff: 2.0,
      exponent: 1.0,
    };
    let mass = law.mass(1.0);
    // 2.0 * 1.0^2 / 2.0 = 1.0
    assert_abs_diff_eq!(1.0, mass, epsilon = 1e-14);
  }

  #[test]
  fn test_approach_law_mass_b2() {
    // p(t) = C * t^2, mass = C * t^3 / 3
    let law = HardApproachLaw {
      t_hard: 0.0,
      coeff: 6.0,
      exponent: 2.0,
    };
    let mass = law.mass(1.0);
    // 6.0 * 1.0^3 / 3.0 = 2.0
    assert_abs_diff_eq!(2.0, mass, epsilon = 1e-14);
  }

  // --- HardApproachLaw::compose_multiply ---

  #[test]
  fn test_approach_law_compose_multiply_exponents_add() {
    let a = HardApproachLaw {
      t_hard: 0.0,
      coeff: 2.0,
      exponent: 1.0,
    };
    let b = HardApproachLaw {
      t_hard: 0.0,
      coeff: 3.0,
      exponent: 2.0,
    };
    let result = a.compose_multiply(&b);
    assert_abs_diff_eq!(3.0, result.exponent, epsilon = 1e-14);
    assert_abs_diff_eq!(6.0, result.coeff, epsilon = 1e-14);
    assert_ulps_eq!(0.0, result.t_hard);
  }

  #[test]
  fn test_approach_law_compose_multiply_b0_b0() {
    let a = HardApproachLaw {
      t_hard: 0.0,
      coeff: 2.0,
      exponent: 0.0,
    };
    let b = HardApproachLaw {
      t_hard: 0.0,
      coeff: 5.0,
      exponent: 0.0,
    };
    let result = a.compose_multiply(&b);
    assert_abs_diff_eq!(0.0, result.exponent, epsilon = 1e-14);
    assert_abs_diff_eq!(10.0, result.coeff, epsilon = 1e-14);
  }

  // --- HardApproachLaw::negate_arg ---

  #[test]
  fn test_approach_law_negate_arg() {
    let law = HardApproachLaw {
      t_hard: 5.0,
      coeff: 2.0,
      exponent: 1.0,
    };
    let negated = law.negate_arg();
    assert_ulps_eq!(-5.0, negated.t_hard);
    assert_abs_diff_eq!(2.0, negated.coeff, epsilon = 1e-14);
    assert_abs_diff_eq!(1.0, negated.exponent, epsilon = 1e-14);
  }

  // --- GridFn integration: approach law in extrapolation ---

  #[test]
  fn test_gridfn_approach_law_left_interpolation() -> Result<(), Report> {
    // Grid [1.0, 5.0], hard left boundary at t=0 with approach law p(t) = 2.0 * t^1.0
    let grid = GridFn::from_range_values((1.0, 5.0), ndarray::array![2.0, 4.0, 6.0, 8.0, 10.0])?.with_left_extrap(
      BoundaryBehavior::Hard(Some(HardApproachLaw {
        t_hard: 0.0,
        coeff: 2.0,
        exponent: 1.0,
      })),
    );

    // Between hard boundary and grid: use approach law
    let val = grid.interp(0.5)?;
    assert_abs_diff_eq!(1.0, val, epsilon = 1e-14); // 2.0 * 0.5^1.0

    // Beyond hard boundary: zero
    let val = grid.interp(-0.5)?;
    assert_abs_diff_eq!(0.0, val, epsilon = 1e-14);

    // On grid: normal interpolation
    let val = grid.interp(1.0)?;
    assert_abs_diff_eq!(2.0, val, epsilon = 1e-14);

    Ok(())
  }

  #[test]
  fn test_gridfn_approach_law_right_interpolation() -> Result<(), Report> {
    // Grid [0.0, 4.0], hard right boundary at t=5.0 with approach law p(t) = 3.0 * (5-t)^2
    let grid = GridFn::from_range_values((0.0, 4.0), ndarray::array![75.0, 48.0, 27.0, 12.0, 3.0])?.with_right_extrap(
      BoundaryBehavior::Hard(Some(HardApproachLaw {
        t_hard: 5.0,
        coeff: 3.0,
        exponent: 2.0,
      })),
    );

    // Between grid and hard boundary: use approach law
    let val = grid.interp(4.5)?;
    assert_abs_diff_eq!(0.75, val, epsilon = 1e-14); // 3.0 * 0.5^2

    // Beyond hard boundary: zero
    let val = grid.interp(5.5)?;
    assert_abs_diff_eq!(0.0, val, epsilon = 1e-14);

    Ok(())
  }

  #[test]
  fn test_gridfn_approach_law_b0_step_at_boundary() -> Result<(), Report> {
    // b=0: density is constant (maximal) at boundary
    let grid = GridFn::from_range_values((0.1, 1.0), ndarray::array![5.0, 5.0, 5.0, 5.0, 5.0])?.with_left_extrap(
      BoundaryBehavior::Hard(Some(HardApproachLaw {
        t_hard: 0.0,
        coeff: 5.0,
        exponent: 0.0,
      })),
    );

    // Between hard boundary and grid: constant value
    let val = grid.interp(0.05)?;
    assert_abs_diff_eq!(5.0, val, epsilon = 1e-14);

    // At hard boundary itself: b=0 returns coeff
    let val = grid.interp(0.0)?;
    assert_abs_diff_eq!(5.0, val, epsilon = 1e-14);

    // Beyond: zero
    let val = grid.interp(-0.1)?;
    assert_abs_diff_eq!(0.0, val, epsilon = 1e-14);

    Ok(())
  }

  #[test]
  fn test_gridfn_scale_y_preserves_approach_law() -> Result<(), Report> {
    let law = HardApproachLaw {
      t_hard: 0.0,
      coeff: 2.0,
      exponent: 1.5,
    };
    let grid = GridFn::from_range_values((1.0, 3.0), ndarray::array![1.0, 2.0, 3.0])?
      .with_left_extrap(BoundaryBehavior::Hard(Some(law)));

    let scaled = grid.scale_y(3.0);
    let scaled_law = scaled
      .left_extrap()
      .approach_law()
      .expect("approach law should be preserved");
    assert_abs_diff_eq!(6.0, scaled_law.coeff, epsilon = 1e-14); // 2.0 * 3.0
    assert_abs_diff_eq!(1.5, scaled_law.exponent, epsilon = 1e-14);
    assert_ulps_eq!(0.0, scaled_law.t_hard);
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_clears_approach_law() -> Result<(), Report> {
    let law = HardApproachLaw {
      t_hard: 0.0,
      coeff: 2.0,
      exponent: 1.0,
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
      coeff: 2.0,
      exponent: 1.0,
    };
    let right_law = HardApproachLaw {
      t_hard: 10.0,
      coeff: 3.0,
      exponent: 2.0,
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

    // Right law (t_hard=10) becomes left law (t_hard=-10)
    assert_ulps_eq!(-10.0, new_left.t_hard);
    assert_abs_diff_eq!(3.0, new_left.coeff, epsilon = 1e-14);
    assert_abs_diff_eq!(2.0, new_left.exponent, epsilon = 1e-14);

    // Left law (t_hard=0) becomes right law (t_hard=0)
    assert_ulps_eq!(0.0, new_right.t_hard);
    assert_abs_diff_eq!(2.0, new_right.coeff, epsilon = 1e-14);
    assert_abs_diff_eq!(1.0, new_right.exponent, epsilon = 1e-14);

    Ok(())
  }

  #[test]
  fn test_gridfn_resample_preserves_approach_law() -> Result<(), Report> {
    let law = HardApproachLaw {
      t_hard: 0.0,
      coeff: 2.0,
      exponent: 1.0,
    };
    let grid = GridFn::from_range_values((1.0, 5.0), ndarray::array![2.0, 4.0, 6.0, 8.0, 10.0])?
      .with_left_extrap(BoundaryBehavior::Hard(Some(law)));

    let resampled = grid.resample_range_dx((1.0, 5.0), 0.5)?;
    let resampled_law = resampled
      .left_extrap()
      .approach_law()
      .expect("approach law should be preserved");
    assert_abs_diff_eq!(2.0, resampled_law.coeff, epsilon = 1e-14);
    assert_abs_diff_eq!(1.0, resampled_law.exponent, epsilon = 1e-14);
    Ok(())
  }
}
