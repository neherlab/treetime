#[cfg(test)]
mod tests {
  use crate::split_fraction::SplitFraction;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::parent_end( 0.0)]
  #[case::middle(     0.5)]
  #[case::child_end(  1.0)]
  #[trace]
  fn test_split_fraction_accepts_values_in_the_unit_interval(#[case] value: f64) {
    let actual = SplitFraction::new(value).unwrap().value();

    assert_eq!(value.to_bits(), actual.to_bits());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::negative(          -0.1,              "A split fraction must be a finite number in [0, 1], got -0.1. This is an internal error. Please report it to developers.")]
  #[case::above_one(         1.1,               "A split fraction must be a finite number in [0, 1], got 1.1. This is an internal error. Please report it to developers.")]
  #[case::nan(               f64::NAN,          "A split fraction must be a finite number in [0, 1], got NaN. This is an internal error. Please report it to developers.")]
  #[case::positive_infinity( f64::INFINITY,     "A split fraction must be a finite number in [0, 1], got inf. This is an internal error. Please report it to developers.")]
  #[case::negative_infinity( f64::NEG_INFINITY, "A split fraction must be a finite number in [0, 1], got -inf. This is an internal error. Please report it to developers.")]
  #[trace]
  fn test_split_fraction_rejects_values_outside_the_unit_interval(#[case] value: f64, #[case] expected: &str) {
    let result = SplitFraction::new(value);

    assert_error!(result, expected);
  }
}
