#[cfg(test)]
mod tests {
  use crate::timetree::inference::tail_fit::fit_message_soft_tail;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_grid::{BoundaryBehavior, Side};
  use treetime_utils::assert_error;

  // A message whose stored neg-log ordinates rise by 1 per unit time decays to the right (higher
  // neg-log = lower probability). The right-side fit recovers that slope exactly: least squares over
  // an exact line is exact. Oracle: the constructed slope.
  #[test]
  fn test_tail_fit_right_side_fits_positive_slope() -> Result<(), Report> {
    let message: Distribution<NegLog> = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![0.0, 1.0, 2.0, 3.0, 4.0])?;
    let BoundaryBehavior::Linear(law) = fit_message_soft_tail(&message, Side::Right)? else {
      panic!("expected a fitted Linear soft tail");
    };
    assert_ulps_eq!(1.0, law.slope, max_ulps = 10);
    Ok(())
  }

  // The mirror case: ordinates falling by 1 per unit time decay to the left, so the left-side fit
  // recovers a slope of -1. Oracle: the constructed slope.
  #[test]
  fn test_tail_fit_left_side_fits_negative_slope() -> Result<(), Report> {
    let message: Distribution<NegLog> = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![4.0, 3.0, 2.0, 1.0, 0.0])?;
    let BoundaryBehavior::Linear(law) = fit_message_soft_tail(&message, Side::Left)? else {
      panic!("expected a fitted Linear soft tail");
    };
    assert_ulps_eq!(-1.0, law.slope, max_ulps = 10);
    Ok(())
  }

  // The property the whole change exists for: a fitted decaying tail is integrable. A flat Constant
  // tail has infinite mass; the Linear tail's half-line integral exp(-y_edge)/|slope| is finite.
  #[test]
  fn test_tail_fit_fitted_tail_has_finite_mass() -> Result<(), Report> {
    let message: Distribution<NegLog> = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![0.0, 1.0, 2.0, 3.0, 4.0])?;
    let BoundaryBehavior::Linear(law) = fit_message_soft_tail(&message, Side::Right)? else {
      panic!("expected a fitted Linear soft tail");
    };
    let mass = law.mass(4.0);
    assert!(mass.is_finite(), "fitted soft-tail mass must be finite, got {mass}");
    assert!(mass > 0.0, "fitted soft-tail mass must be positive, got {mass}");
    Ok(())
  }

  // Point, Range, and Empty messages carry no grid tail; the helper returns an inert Constant
  // because with_*_extrap is a no-op on non-Function distributions.
  #[test]
  fn test_tail_fit_point_message_is_inert_constant() -> Result<(), Report> {
    let message: Distribution<NegLog> = Distribution::point(2015.0, 0.0);
    assert_eq!(BoundaryBehavior::Constant, fit_message_soft_tail(&message, Side::Left)?);
    Ok(())
  }

  #[test]
  fn test_tail_fit_range_message_is_inert_constant() -> Result<(), Report> {
    let message: Distribution<NegLog> = Distribution::range((2014.0, 2015.0), 0.0);
    assert_eq!(BoundaryBehavior::Constant, fit_message_soft_tail(&message, Side::Right)?);
    Ok(())
  }

  #[test]
  fn test_tail_fit_empty_message_is_inert_constant() -> Result<(), Report> {
    let message: Distribution<NegLog> = Distribution::empty();
    assert_eq!(BoundaryBehavior::Constant, fit_message_soft_tail(&message, Side::Left)?);
    Ok(())
  }

  // A grid with fewer than two finite points in the fit window cannot yield a slope. The helper
  // errors rather than fabricating a flat fallback, which would silently reintroduce the
  // non-integrable Constant behavior. Here only the right edge is finite.
  #[test]
  fn test_tail_fit_degenerate_grid_errors() {
    let message: Distribution<NegLog> =
      Distribution::function(array![0.0, 1.0], array![f64::INFINITY, 0.0]).expect("two-point function");
    assert_error!(
      fit_message_soft_tail(&message, Side::Right),
      "Timetree message cannot fit a soft tail on the Right side: the convolved grid is degenerate. This is an internal error. Please report it to developers."
    );
  }
}
