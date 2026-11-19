#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionNegLog;
  use crate::distribution_ops::negate::distribution_negation;
  use eyre::Report;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::empty(          DistributionNegLog::empty(),               DistributionNegLog::empty())]
  #[case::point(          DistributionNegLog::point(2.0, 3.0),       DistributionNegLog::point(-2.0, 3.0))]
  #[case::range_reversed( DistributionNegLog::range((1.0, 4.0), 2.0), DistributionNegLog::range((-4.0, -1.0), 2.0))]
  #[case::function_mirrored(
    DistributionNegLog::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap(),
    DistributionNegLog::function(array![-2.0, -1.0, 0.0], array![3.0, 2.0, 1.0]).unwrap(),
  )]
  #[trace]
  fn test_negation_mirrors_the_distribution_about_zero(
    #[case] dist: DistributionNegLog,
    #[case] expected: DistributionNegLog,
  ) -> Result<(), Report> {
    let actual = distribution_negation(&dist)?;
    assert_eq!(expected, actual);
    Ok(())
  }
}
