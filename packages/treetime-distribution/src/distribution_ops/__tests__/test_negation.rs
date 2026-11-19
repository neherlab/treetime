#[cfg(test)]
mod tests {
  use crate::Distribution;
  use crate::distribution_ops::negate::distribution_negation;
  use eyre::Report;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::empty(          Distribution::empty(),               Distribution::empty())]
  #[case::point(          Distribution::point(2.0, 3.0),       Distribution::point(-2.0, 3.0))]
  #[case::range_reversed( Distribution::range((1.0, 4.0), 2.0), Distribution::range((-4.0, -1.0), 2.0))]
  #[case::function_mirrored(
    Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap(),
    Distribution::function(array![-2.0, -1.0, 0.0], array![3.0, 2.0, 1.0]).unwrap(),
  )]
  #[trace]
  fn test_negation_mirrors_the_distribution_about_zero(
    #[case] dist: Distribution,
    #[case] expected: Distribution,
  ) -> Result<(), Report> {
    let actual = distribution_negation(&dist)?;
    assert_eq!(expected, actual);
    Ok(())
  }
}
