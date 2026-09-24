#[cfg(test)]
mod tests {
  use crate::fmt::float::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  #[case::decimal_truncated((1.23456, 3), "1.23")]
  #[case::large_integer((123.456, 3), "123")]
  #[case::small_decimal((0.00123, 3), "0.00123")]
  #[case::round_thousand((1000.0,  3), "1000")]
  #[trace]
  fn test_float_to_significant_digits(#[case] (input, digits): (f64, u8), #[case] expected: &str) {
    assert_eq!(float_to_significant_digits(input, digits), expected);
  }

  #[rstest]
  #[case::truncate_decimals((1.23456, 2), "1.23")]
  #[case::round_up((123.456, 2), "123.46")]
  #[case::trim_trailing_zeros((1.0,     2), "1")]
  #[case::zero_decimals((1.0,     0), "1")]
  #[trace]
  fn test_float_to_decimal_digits(#[case] (input, digits): (f64, i8), #[case] expected: &str) {
    assert_eq!(float_to_decimal_digits(input, digits), expected);
  }

  #[rstest]
  #[case::sig_only((1.23456, Some(3), None   ), "1.23")]
  #[case::dec_only((1.23456, None,    Some(2)), "1.23")]
  #[case::both_constraints((1.23456, Some(3), Some(2)), "1.23")]
  #[case::default_behavior((1.23456, None,    None   ), "1.23")]
  #[trace]
  fn test_float_to_digits(
    #[case] (input, sig_digits, dec_digits): (f64, Option<u8>, Option<i8>),
    #[case] expected: &str,
  ) {
    assert_eq!(float_to_digits(input, sig_digits, dec_digits), expected);
  }

  #[rstest]
  #[case::f64_decimal((1.23456_f64, 3), "1.23")]
  #[case::f32_integer((123.456_f32, 3), "123")]
  #[case::f64_small((0.00123_f64, 3), "0.00123")]
  #[trace]
  fn test_trait_to_significant_digits(#[case] (input, digits): (impl FloatFormatExt, u8), #[case] expected: &str) {
    assert_eq!(input.to_significant_digits(digits), expected);
  }

  #[rstest]
  #[case::f64_truncate((1.23456_f64, 2), "1.23")]
  #[case::f32_round_up((123.456_f32, 2), "123.46")]
  #[case::f64_trailing_zeros((1.0_f64,     2), "1")]
  #[trace]
  fn test_trait_to_decimal_digits(#[case] (input, digits): (impl FloatFormatExt, i8), #[case] expected: &str) {
    assert_eq!(input.to_decimal_digits(digits), expected);
  }

  #[rstest]
  #[case::f64_sig_only((1.23456_f64, Some(3), None   ), "1.23")]
  #[case::f32_dec_only((1.23456_f32, None,    Some(2)), "1.23")]
  #[case::f64_both((1.23456_f64, Some(3), Some(2)), "1.23")]
  #[case::f32_default((1.23456_f32, None,    None   ), "1.23")]
  #[trace]
  fn test_trait_to_digits(
    #[case] (input, sig_digits, dec_digits): (impl FloatFormatExt, Option<u8>, Option<i8>),
    #[case] expected: &str,
  ) {
    assert_eq!(input.to_digits(sig_digits, dec_digits), expected);
  }

  #[rstest]
  #[case::trailing_zero(("1.230", "1.23"))]
  #[case::all_zeros(("1.000", "1"))]
  #[case::single_zero(("1.0",   "1"))]
  #[case::no_decimal(("123",   "123"))]
  #[case::leading_zero(("0.100", "0.1"))]
  #[case::sci_negative_exp(("8.7110e-10", "8.7110e-10"))]
  #[case::sci_positive_exp(("1.2300e+10", "1.2300e+10"))]
  #[case::sci_no_trim(("1.000e-10", "1.000e-10"))]
  #[trace]
  fn test_trim_trailing_zeros(#[case] (input, expected): (&str, &str)) {
    assert_eq!(trim_trailing_zeros(input), expected);
  }

  #[rstest]
  #[case::zero((0.0_f64,           3), "0")]
  #[case::negative((-1.23456_f64,      3), "-1.23")]
  #[case::positive_infinity((f64::INFINITY,     3), "inf")]
  #[case::negative_infinity((f64::NEG_INFINITY, 3), "-inf")]
  #[trace]
  fn test_edge_cases(#[case] (input, digits): (f64, u8), #[case] expected: &str) {
    assert_eq!(input.to_significant_digits(digits), expected);
  }
}
