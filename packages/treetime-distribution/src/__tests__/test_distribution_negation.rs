#[cfg(test)]
mod tests {
  use crate::DistributionPlain as Distribution;
  use crate::distribution_negation::{distribution_negation, distribution_negation_inplace};
  use eyre::Report;
  use ndarray::array;

  #[test]
  fn test_negate_empty() -> Result<(), Report> {
    let dist: Distribution = Distribution::empty();
    let actual: Distribution = distribution_negation(&dist)?;
    let expected: Distribution = Distribution::empty();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_point() -> Result<(), Report> {
    let dist: Distribution = Distribution::point(2.0, 3.0);
    let actual: Distribution = distribution_negation(&dist)?;
    let expected: Distribution = Distribution::point(-2.0, 3.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_point_zero() -> Result<(), Report> {
    let dist: Distribution = Distribution::point(0.0, 5.0);
    let actual: Distribution = distribution_negation(&dist)?;
    let expected: Distribution = Distribution::point(0.0, 5.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_range() -> Result<(), Report> {
    let dist: Distribution = Distribution::range((1.0, 4.0), 2.0);
    let actual: Distribution = distribution_negation(&dist)?;
    let expected: Distribution = Distribution::range((-4.0, -1.0), 2.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_range_symmetric() -> Result<(), Report> {
    let dist: Distribution = Distribution::range((-3.0, 3.0), 1.0);
    let actual: Distribution = distribution_negation(&dist)?;
    let expected: Distribution = Distribution::range((-3.0, 3.0), 1.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_function() -> Result<(), Report> {
    let t = array![0.0, 1.0, 2.0];
    let y = array![1.0, 2.0, 3.0];
    let dist: Distribution = Distribution::function(t, y)?;

    let actual: Distribution = distribution_negation(&dist)?;

    let expected_t = array![-2.0, -1.0, 0.0];
    let expected_y = array![3.0, 2.0, 1.0];
    let expected: Distribution = Distribution::function(expected_t, expected_y)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_inplace_point() -> Result<(), Report> {
    let mut actual: Distribution = Distribution::point(2.0, 3.0);
    distribution_negation_inplace(&mut actual)?;
    let expected: Distribution = Distribution::point(-2.0, 3.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_inplace_range() -> Result<(), Report> {
    let mut actual: Distribution = Distribution::range((1.0, 4.0), 2.0);
    distribution_negation_inplace(&mut actual)?;
    let expected: Distribution = Distribution::range((-4.0, -1.0), 2.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_negate_inplace_function() -> Result<(), Report> {
    let t = array![0.0, 1.0, 2.0];
    let y = array![1.0, 2.0, 3.0];
    let mut actual: Distribution = Distribution::function(t, y)?;

    distribution_negation_inplace(&mut actual)?;

    let expected_t = array![-2.0, -1.0, 0.0];
    let expected_y = array![3.0, 2.0, 1.0];
    let expected: Distribution = Distribution::function(expected_t, expected_y)?;
    assert_eq!(expected, actual);
    Ok(())
  }
}
