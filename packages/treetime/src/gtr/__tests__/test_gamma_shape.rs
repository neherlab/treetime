#[cfg(test)]
mod tests {
  use crate::gtr::gamma::{GAMMA_SHAPE_MAX, GAMMA_SHAPE_MIN};
  use crate::gtr::gamma_shape::{ShapeEstimate, optimize_gamma_shape};
  use crate::make_error;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_primitives::LogLh;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::start_at_optimum(        0.7,  0.7)]
  #[case::start_above_optimum(     0.7,  50.0)]
  #[case::start_below_optimum(     0.7,  GAMMA_SHAPE_MIN)]
  #[case::start_at_upper_bound(    3.0,  GAMMA_SHAPE_MAX)]
  #[case::optimum_near_lower_bound(0.021, 1.0)]
  #[trace]
  fn test_gamma_shape_finds_the_interior_maximum(#[case] optimum: f64, #[case] start: f64) -> Result<(), Report> {
    let ShapeEstimate { shape, .. } = optimize_gamma_shape(start, |shape| Ok(log_lh_peaked_at(optimum, shape)))?;
    pretty_assert_abs_diff_eq!(optimum.ln(), shape.ln(), epsilon = 1e-4);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::below_range_from_inside( 0.001, 1.0,             GAMMA_SHAPE_MIN)]
  #[case::below_range_from_bound(  0.001, GAMMA_SHAPE_MIN, GAMMA_SHAPE_MIN)]
  #[case::above_range_from_inside( 500.0, 1.0,             GAMMA_SHAPE_MAX)]
  #[case::above_range_from_bound(  500.0, GAMMA_SHAPE_MAX, GAMMA_SHAPE_MAX)]
  #[trace]
  fn test_gamma_shape_stops_at_the_bound_beyond_which_the_maximum_lies(
    #[case] optimum: f64,
    #[case] start: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let ShapeEstimate { shape, .. } = optimize_gamma_shape(start, |shape| Ok(log_lh_peaked_at(optimum, shape)))?;
    pretty_assert_ulps_eq!(expected, shape, max_ulps = 0);
    Ok(())
  }

  #[test]
  fn test_gamma_shape_at_the_bound_needs_two_evaluations_from_the_bound() -> Result<(), Report> {
    let estimate = optimize_gamma_shape(GAMMA_SHAPE_MIN, |shape| Ok(log_lh_peaked_at(0.001, shape)))?;
    pretty_assert_ulps_eq!(GAMMA_SHAPE_MIN, estimate.shape, max_ulps = 0);
    assert_eq!(2, estimate.evaluations);
    Ok(())
  }

  #[test]
  fn test_gamma_shape_reports_a_failed_evaluation() {
    let result = optimize_gamma_shape(1.0, |shape| {
      if shape > 1.0 {
        return make_error!("likelihood failed at {shape}");
      }
      Ok(log_lh_peaked_at(0.5, shape))
    });
    assert_error!(result, "likelihood failed at 1.1051709180756477");
  }

  fn log_lh_peaked_at(optimum: f64, shape: f64) -> LogLh {
    LogLh::new(-500.0 * (shape.ln() - optimum.ln()).powi(2))
  }
}
