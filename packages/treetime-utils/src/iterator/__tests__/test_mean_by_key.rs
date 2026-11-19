#[cfg(test)]
mod tests {
  use crate::iterator::mean_by_key::MeanByKey;
  use crate::pretty_assert_ulps_eq;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::empty(    vec![],                  0)]
  #[case::single(   vec![7],                 14)]
  #[case::several(  vec![1, 2, 3, 4, 5],     6)]
  #[case::negative( vec![-5, -2, -1, 0, 3],  -2)]
  #[trace]
  fn test_mean_by_key_returns_mean_of_keys(#[case] data: Vec<i32>, #[case] expected: i32) {
    let actual: i32 = data.iter().mean_by_key(|x| *x * 2);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_mean_by_key_returns_fractional_mean_of_float_keys() {
    let data = [(1.0, "a"), (2.0, "b"), (3.0, "c"), (4.0, "d")];
    let actual: f64 = data.iter().mean_by_key(|(x, _)| *x);
    pretty_assert_ulps_eq!(2.5, actual);
  }
}
