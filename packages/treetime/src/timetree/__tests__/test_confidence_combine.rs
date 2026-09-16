#[cfg(test)]
mod tests {
  use crate::timetree::confidence::combine_confidence;
  use approx::assert_relative_eq;

  #[test]
  fn test_combine_confidence_no_contributions() {
    let result = combine_confidence(10.0, (5.0, 15.0), None, None);
    assert_relative_eq!(result.0, 5.0);
    assert_relative_eq!(result.1, 15.0);
  }

  #[test]
  fn test_combine_confidence_single_contribution() {
    let result = combine_confidence(10.0, (0.0, 20.0), Some((8.0, 12.0)), None);
    assert_relative_eq!(result.0, 8.0);
    assert_relative_eq!(result.1, 12.0);

    let result = combine_confidence(10.0, (0.0, 20.0), None, Some((7.0, 13.0)));
    assert_relative_eq!(result.0, 7.0);
    assert_relative_eq!(result.1, 13.0);
  }

  #[test]
  fn test_combine_confidence_quadrature() {
    // c1: (8, 12) -> deviations of 2 from center 10
    // c2: (7, 13) -> deviations of 3 from center 10
    // Combined: sqrt(2^2 + 3^2) = sqrt(13) = 3.606
    let result = combine_confidence(10.0, (0.0, 20.0), Some((8.0, 12.0)), Some((7.0, 13.0)));
    let expected_dev = (4.0_f64 + 9.0).sqrt();
    assert_relative_eq!(result.0, 10.0 - expected_dev, epsilon = 1e-10);
    assert_relative_eq!(result.1, 10.0 + expected_dev, epsilon = 1e-10);
  }

  #[test]
  fn test_combine_confidence_clipped_to_limits() {
    // Large contributions that exceed limits
    let result = combine_confidence(10.0, (5.0, 15.0), Some((0.0, 20.0)), Some((0.0, 20.0)));
    // Quadrature would give sqrt(100 + 100) = 14.14 deviation
    // But limits clip to (5, 15)
    assert_relative_eq!(result.0, 5.0);
    assert_relative_eq!(result.1, 15.0);
  }
}
