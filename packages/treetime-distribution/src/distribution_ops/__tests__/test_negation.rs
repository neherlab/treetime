#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionPlain;
  use crate::distribution_ops::negate::distribution_negation;
  use eyre::Report;
  use ndarray::array;

  #[test]
  fn test_negate_empty() -> Result<(), Report> {
    let dist: DistributionPlain = DistributionPlain::empty();
    let actual: DistributionPlain = distribution_negation(&dist)?;
    let expected: DistributionPlain = DistributionPlain::empty();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_point() -> Result<(), Report> {
    let dist: DistributionPlain = DistributionPlain::point(2.0, 3.0);
    let actual: DistributionPlain = distribution_negation(&dist)?;
    let expected: DistributionPlain = DistributionPlain::point(-2.0, 3.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_point_zero() -> Result<(), Report> {
    let dist: DistributionPlain = DistributionPlain::point(0.0, 5.0);
    let actual: DistributionPlain = distribution_negation(&dist)?;
    let expected: DistributionPlain = DistributionPlain::point(0.0, 5.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_range() -> Result<(), Report> {
    let dist: DistributionPlain = DistributionPlain::range((1.0, 4.0), 2.0);
    let actual: DistributionPlain = distribution_negation(&dist)?;
    let expected: DistributionPlain = DistributionPlain::range((-4.0, -1.0), 2.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_range_symmetric() -> Result<(), Report> {
    let dist: DistributionPlain = DistributionPlain::range((-3.0, 3.0), 1.0);
    let actual: DistributionPlain = distribution_negation(&dist)?;
    let expected: DistributionPlain = DistributionPlain::range((-3.0, 3.0), 1.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_function() -> Result<(), Report> {
    let t = array![0.0, 1.0, 2.0];
    let y = array![1.0, 2.0, 3.0];
    let dist: DistributionPlain = DistributionPlain::function(t, y)?;

    let actual: DistributionPlain = distribution_negation(&dist)?;

    let expected_t = array![-2.0, -1.0, 0.0];
    let expected_y = array![3.0, 2.0, 1.0];
    let expected: DistributionPlain = DistributionPlain::function(expected_t, expected_y)?;
    assert_eq!(expected, actual);
    Ok(())
  }
}
