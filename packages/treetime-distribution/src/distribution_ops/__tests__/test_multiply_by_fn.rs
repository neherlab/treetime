#[cfg(test)]
mod tests {
  use self::helpers::{flat, linear_with_soft_side, sides};
  use crate::__tests__::aliases::DistributionNegLog;
  use crate::DistributionFunction;
  use crate::distribution_ops::multiply_by_fn::distribution_multiply_by_fn;
  use eyre::Report;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_grid::{BoundaryBehavior, Side, SoftTailLaw};
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  #[test]
  fn test_multiply_by_fn_constant_is_uniform_shift() {
    let neglog = DistributionNegLog::function(array![0.0, 1.0, 2.0], array![4.0, 0.0, 3.0]).unwrap();
    let actual = distribution_multiply_by_fn(&neglog, |_| Ok(10.0), flat).unwrap();
    pretty_assert_abs_diff_eq!(array![14.0, 10.0, 13.0], actual.y().unwrap(), epsilon = 1e-15);
  }

  #[test]
  fn test_multiply_by_fn_empty_passthrough() {
    let actual = distribution_multiply_by_fn(&DistributionNegLog::Empty, |_| Ok(1.0), flat).unwrap();
    assert_eq!(DistributionNegLog::Empty, actual);
  }

  #[test]
  fn test_multiply_by_fn_rejects_formula() {
    let formula = DistributionNegLog::Formula(crate::DistributionFormula::new(|_| Ok(1.0), 0.0, 10.0));
    assert_error!(
      distribution_multiply_by_fn(&formula, |_| Ok(1.0), flat),
      "distribution_multiply_by_fn requires a concrete Point, Range, or Function distribution"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::left( Side::Left,  -2.0,  0.5, (-1.0, 1.5))]
  #[case::right(Side::Right,  2.0, -0.5, ( 5.0, 7.5))]
  #[trace]
  fn test_multiply_by_fn_adds_the_weight_slope_beyond_the_soft_edge(
    #[case] side: Side,
    #[case] slope: f64,
    #[case] weight_slope: f64,
    #[case] (t, expected): (f64, f64),
  ) -> Result<(), Report> {
    let input = linear_with_soft_side(slope, side)?;
    let at_the_edge = |at: Side, edge: f64| {
      pretty_assert_ulps_eq!(if at == Side::Left { 0.0 } else { 4.0 }, edge, max_ulps = 0);
      weight_slope
    };

    let actual = distribution_multiply_by_fn(&input, |ti: f64| Ok(weight_slope * ti), at_the_edge)?;

    pretty_assert_abs_diff_eq!(expected, actual.eval(t)?, epsilon = 1e-14);
    let DistributionNegLog::Function(rf) = actual else {
      panic!("expected a Function result, got {actual:?}");
    };
    let soft = BoundaryBehavior::Linear(SoftTailLaw { slope: slope + weight_slope });
    let expected_laws = match side {
      Side::Left => (soft, BoundaryBehavior::Hard),
      Side::Right => (BoundaryBehavior::Hard, soft),
    };
    assert_eq!(expected_laws, (rf.left_extrap(), rf.right_extrap()));
    Ok(())
  }

  #[test]
  fn test_multiply_by_fn_step_one_cell_inside_the_edge_keeps_a_decaying_tail() -> Result<(), Report> {
    let input = linear_with_soft_side(-2.0, Side::Left)?;
    let weight = |ti: f64| if ti < 0.2 { 0.0 } else { 3.0 };

    let actual = distribution_multiply_by_fn(&input, |ti: f64| Ok(weight(ti)), flat)?;

    let t = Array1::linspace(0.0, 4.0, 21);
    let p = t.mapv(|ti| (2.0 * ti - weight(ti)).exp());
    let grid_mass = 0.2 * (p.sum() - f64::midpoint(p[0], p[20]));
    let tail_mass = 1.0 / 2.0;
    pretty_assert_abs_diff_eq!((grid_mass + tail_mass).ln(), actual.log_mass()?, epsilon = 1e-14);
    let DistributionNegLog::Function(rf) = actual else {
      panic!("expected a Function result, got {actual:?}");
    };
    assert_eq!(BoundaryBehavior::Linear(SoftTailLaw { slope: -2.0 }), rf.left_extrap());
    Ok(())
  }

  #[test]
  fn test_multiply_by_fn_tail_law_takes_the_weight_slope_at_the_edge_and_omits_a_step_farther_out() -> Result<(), Report>
  {
    let incoming_slope = -1.0;
    let input = linear_with_soft_side(incoming_slope, Side::Left)?;
    let weight = |ti: f64| if ti < -0.5 { -10.0 } else { 0.0 };
    let slope_at_edge = 0.0;

    let actual = distribution_multiply_by_fn(&input, |ti: f64| Ok(weight(ti)), sides(slope_at_edge, 0.0))?;

    let composed_slope = incoming_slope + slope_at_edge;
    let (edge, t) = (0.0, -1.0);
    let ordinate_at_edge = incoming_slope * edge + weight(edge);
    let expected = ordinate_at_edge + composed_slope * (t - edge);
    pretty_assert_abs_diff_eq!(expected, actual.eval(t)?, epsilon = 1e-15);
    let DistributionNegLog::Function(rf) = actual else {
      panic!("expected a Function result, got {actual:?}");
    };
    assert_eq!(
      BoundaryBehavior::Linear(SoftTailLaw { slope: composed_slope }),
      rf.left_extrap()
    );
    Ok(())
  }

  #[test]
  fn test_multiply_by_fn_log_mass_counts_the_composed_tail() -> Result<(), Report> {
    let incoming_slope = -2.0;
    let input = linear_with_soft_side(incoming_slope, Side::Left)?;
    let (slope_inside, slope_outside) = (2.0, 1.5);
    let weight = move |ti: f64| {
      if ti < 0.0 {
        slope_outside * ti
      } else {
        slope_inside * ti
      }
    };
    let actual = distribution_multiply_by_fn(&input, |ti: f64| Ok(weight(ti)), sides(slope_outside, slope_inside))?;

    let grid_mass: f64 = 4.0;
    let tail_mass = 1.0 / -(incoming_slope + slope_outside);
    pretty_assert_abs_diff_eq!((grid_mass + tail_mass).ln(), actual.log_mass()?, epsilon = 1e-15);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::left( Side::Left,  -2.0,  3.0, "The weight does not let the soft Left tail decay: incoming slope -2, weight slope 3 beyond the edge at 0. This is an internal error. Please report it to developers.")]
  #[case::right(Side::Right,  2.0, -2.0, "The weight does not let the soft Right tail decay: incoming slope 2, weight slope -2 beyond the edge at 4. This is an internal error. Please report it to developers.")]
  #[trace]
  fn test_multiply_by_fn_rejects_a_weight_that_stops_the_tail_decay(
    #[case] side: Side,
    #[case] slope: f64,
    #[case] weight_slope: f64,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let input = linear_with_soft_side(slope, side)?;
    assert_error!(
      distribution_multiply_by_fn(&input, |ti: f64| Ok(weight_slope * ti), sides(weight_slope, weight_slope)),
      expected
    );
    Ok(())
  }

  #[test]
  fn test_multiply_by_fn_carries_hard_both_sides() {
    let t = Array1::linspace(0.0, 4.0, 21);
    let y = t.mapv(|ti| -2.0 * ti);
    let input = DistributionFunction::from_arrays(&t, y)
      .unwrap()
      .with_left_extrap(BoundaryBehavior::Hard)
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();
    let actual = distribution_multiply_by_fn(
      &DistributionNegLog::Function(input),
      |ti: f64| Ok(0.5 * ti),
      sides(0.5, 0.5),
    )
    .unwrap();

    let DistributionNegLog::Function(rf) = actual else {
      panic!("expected a Function result");
    };
    assert_eq!(BoundaryBehavior::Hard, rf.left_extrap());
    assert_eq!(BoundaryBehavior::Hard, rf.right_extrap());
  }

  mod helpers {
    use crate::__tests__::aliases::DistributionNegLog;
    use crate::DistributionFunction;
    use eyre::Report;
    use ndarray::Array1;
    use treetime_grid::{BoundaryBehavior, Side, SoftTailLaw};

    pub(super) fn flat(_side: Side, _edge: f64) -> f64 {
      0.0
    }

    pub(super) fn sides(left: f64, right: f64) -> impl Fn(Side, f64) -> f64 {
      move |side, _edge| match side {
        Side::Left => left,
        Side::Right => right,
      }
    }

    pub(super) fn linear_with_soft_side(slope: f64, side: Side) -> Result<DistributionNegLog, Report> {
      let t = Array1::linspace(0.0, 4.0, 21);
      let f = DistributionFunction::from_arrays(&t, t.mapv(|ti| slope * ti))?;
      let soft = BoundaryBehavior::Linear(SoftTailLaw { slope });
      let f = match side {
        Side::Left => f.with_left_extrap(soft)?.with_right_extrap(BoundaryBehavior::Hard)?,
        Side::Right => f.with_left_extrap(BoundaryBehavior::Hard)?.with_right_extrap(soft)?,
      };
      Ok(DistributionNegLog::Function(f))
    }
  }
}
