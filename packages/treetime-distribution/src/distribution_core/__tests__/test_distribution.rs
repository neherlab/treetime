#[cfg(test)]
mod tests {
  use crate::distribution_core::formula::DistributionFormula;
  use crate::Distribution;
  use crate::policy::{NegLog, Plain};
  use ndarray::array;
  use treetime_utils::{assert_error, make_error, pretty_assert_ulps_eq};

  #[test]
  fn test_distribution_time_bounds_empty_is_none() {
    assert_eq!(None, Distribution::<Plain>::Empty.time_bounds());
  }

  #[test]
  fn test_distribution_time_bounds_point_is_degenerate_interval() {
    assert_eq!(Some((1.0, 1.0)), Distribution::<Plain>::point(1.0, 5.0).time_bounds());
  }

  #[test]
  fn test_distribution_time_bounds_range_spans_endpoints() {
    assert_eq!(
      Some((0.0, 2.0)),
      Distribution::<Plain>::range((0.0, 2.0), 3.0).time_bounds()
    );
  }

  #[test]
  fn test_distribution_time_bounds_function_spans_grid() {
    let func = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![1.0, 4.0, 2.0]).unwrap();
    assert_eq!(Some((0.0, 2.0)), func.time_bounds());
  }

  #[test]
  fn test_distribution_y_formula_evaluates_bounds() {
    let formula = Distribution::<NegLog>::Formula(DistributionFormula::new(|t| Ok(2.0 * t), 1.0, 3.0));
    assert_eq!(array![2.0, 6.0], formula.y().unwrap());
  }

  #[test]
  fn test_distribution_y_formula_propagates_evaluation_error() {
    let formula = Distribution::<NegLog>::Formula(DistributionFormula::new(
      |t| make_error!("no value at {t}"),
      1.0,
      3.0,
    ));
    assert_error!(formula.y(), "When evaluating a formula distribution at its bounds [1, 3]: no value at 1");
  }

  #[test]
  fn test_distribution_neglog_normalize_function_shifts_peak_to_zero() {
    let distribution = Distribution::<NegLog>::function(array![0.0, 1.0, 2.0], array![1004.0, 1000.0, 1003.0]).unwrap();
    let normalized = distribution.normalize();
    let expected = array![4.0, 0.0, 3.0];
    let Distribution::Function(f) = normalized else {
      panic!("Expected Function");
    };
    pretty_assert_ulps_eq!(expected, f.y(), max_ulps = 4);
  }

  #[test]
  fn test_distribution_neglog_normalize_is_idempotent() {
    let distribution = Distribution::<NegLog>::function(array![0.0, 1.0, 2.0], array![4.0, 0.0, 3.0]).unwrap();
    let once = distribution.normalize();
    let twice = once.normalize();
    assert_eq!(once, twice);
  }

  #[test]
  fn test_distribution_neglog_normalize_point_maps_to_unit_probability() {
    let normalized = Distribution::<NegLog>::point(2.0, 1000.0).normalize();
    let expected = Distribution::<NegLog>::point(2.0, 0.0);
    assert_eq!(expected, normalized);
  }

  #[test]
  fn test_distribution_neglog_normalize_range_maps_to_unit_probability() {
    let normalized = Distribution::<NegLog>::range((1.0, 3.0), 1000.0).normalize();
    let expected = Distribution::<NegLog>::range((1.0, 3.0), 0.0);
    assert_eq!(expected, normalized);
  }

  #[test]
  fn test_distribution_neglog_normalize_empty_stays_empty() {
    assert_eq!(Distribution::<NegLog>::Empty, Distribution::<NegLog>::Empty.normalize());
  }

  #[test]
  fn test_distribution_neglog_normalize_rejects_nonfinite_minimum() {
    let distribution = Distribution::<NegLog>::function(
      array![0.0, 1.0, 2.0],
      array![f64::INFINITY, f64::INFINITY, f64::INFINITY],
    )
    .unwrap();
    assert_eq!(Distribution::<NegLog>::Empty, distribution.normalize());
  }
}
