#[cfg(test)]
mod tests {
  use crate::iterator::mean_by_key::MeanByKey;
  use crate::{assert_error, pretty_assert_ulps_eq};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::single(   vec![7],                 14)]
  #[case::several(  vec![1, 2, 3, 4, 5],     6)]
  #[case::negative( vec![-5, -2, -1, 0, 3],  -2)]
  #[trace]
  fn test_mean_by_key_returns_mean_of_keys(#[case] data: Vec<i32>, #[case] expected: i32) -> Result<(), Report> {
    let actual: i32 = data.iter().mean_by_key(|x| *x * 2)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_mean_by_key_returns_fractional_mean_of_float_keys() -> Result<(), Report> {
    let data = [(1.0, "a"), (2.0, "b"), (3.0, "c"), (4.0, "d")];
    let actual: f64 = data.iter().mean_by_key(|(x, _)| *x)?;
    pretty_assert_ulps_eq!(2.5, actual);
    Ok(())
  }

  #[test]
  fn test_mean_by_key_rejects_an_empty_sequence() {
    let data: Vec<f64> = vec![];
    let actual = data.iter().mean_by_key(|x| *x);
    assert_error!(actual, "Cannot compute the mean of an empty sequence");
  }

  #[test]
  fn test_mean_by_key_accepts_the_largest_count_that_the_key_type_holds() -> Result<(), Report> {
    let data = [1_i8; 127];
    let actual = data.iter().mean_by_key(|x| *x)?;
    assert_eq!(1_i8, actual);
    Ok(())
  }

  #[test]
  fn test_mean_by_key_rejects_a_count_that_the_key_type_cannot_hold() {
    let data = [0_i8; 128];
    let actual = data.iter().mean_by_key(|x| *x);
    assert_error!(
      actual,
      "Cannot compute the mean of 128 values: the count does not fit in the key type i8"
    );
  }
}
