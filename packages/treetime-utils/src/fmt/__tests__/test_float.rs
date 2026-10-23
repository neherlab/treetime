#[cfg(test)]
mod tests {
  use crate::fmt::float::{
    FloatFormatExt, float_to_digits, float_to_exponential, float_to_power_of_ten, float_to_power_of_ten_with_error,
    float_to_significant_digits, trim_trailing_zeros,
  };
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

  #[rustfmt::skip]
  #[rstest]
  #[case::small_positive(    (0.066_17,   3), "6.617e-02")]
  #[case::negative(          (-70.77,     3), "-7.077e+01")]
  #[case::zero(              (0.0,        3), "0.000e+00")]
  #[case::three_digit_exp(   (1e-100,     3), "1.000e-100")]
  #[case::rounds_mantissa(   (123_456.0,  2), "1.23e+05")]
  #[case::one(               (1.0,        0), "1e+00")]
  #[trace]
  fn test_float_to_exponential(#[case] (input, decimals): (f64, usize), #[case] expected: &str) {
    assert_eq!(float_to_exponential(input, decimals), expected);
  }

  #[rstest]
  #[case::sig_only(            (12.3456, Some(4), None   ), "12.35")]
  #[case::dec_only(            (12.3456, None,    Some(3)), "12.346")]
  #[case::both_decimal_stricter((12.3456, Some(5), Some(1)), "12.3")]
  #[case::both_sig_stricter(   (12.3456, Some(2), Some(3)), "12")]
  #[case::default_behavior(    (12.3456, None,    None   ), "12.3")]
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
  #[case::f64_sig_only(             (12.3456_f64, Some(4), None   ), "12.35")]
  #[case::f32_dec_only(             (12.3456_f32, None,    Some(3)), "12.346")]
  #[case::f64_both_decimal_stricter((12.3456_f64, Some(5), Some(1)), "12.3")]
  #[case::f32_both_sig_stricter(    (12.3456_f32, Some(2), Some(3)), "12")]
  #[case::f32_default(              (12.3456_f32, None,    None   ), "12.3")]
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

  #[rustfmt::skip]
  #[rstest]
  #[case::clock_rate(        (5.072_949e-5,       3), "5.07 × 10⁻⁵")]
  #[case::keeps_zeros(       (3.0e-3,             3), "3.00 × 10⁻³")]
  #[case::rounds_up_a_decade((9.996e-5,           3), "1.00 × 10⁻⁴")]
  #[case::positive_exponent( (6.02e23,            3), "6.02 × 10²³")]
  #[case::negative_minus(    (-2.1e-5,            2), "\u{2212}2.1 × 10⁻⁵")]
  #[case::exponent_zero(     (4.25,               3), "4.25")]
  #[case::single_digit(      (5.0e-5,             1), "5 × 10⁻⁵")]
  #[case::zero(              (0.0,                3), "0")]
  #[case::nan(               (f64::NAN,           3), "NaN")]
  #[case::infinity(          (f64::INFINITY,      3), "∞")]
  #[case::negative_infinity( (f64::NEG_INFINITY,  3), "\u{2212}∞")]
  #[trace]
  fn test_float_to_power_of_ten(#[case] (input, digits): (f64, usize), #[case] expected: &str) {
    assert_eq!(expected, float_to_power_of_ten(input, digits));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::shares_the_exponent((5.072_949e-5, 3.1e-6, 3), "5.07 ± 0.31 × 10⁻⁵")]
  #[case::exponent_zero(      (4.25,         0.126,  3), "4.25 ± 0.13")]
  #[case::zero_value(         (0.0,          3.1e-6, 2), "0 ± 3.1 × 10⁻⁶")]
  #[trace]
  fn test_float_to_power_of_ten_with_error(#[case] (input, error, digits): (f64, f64, usize), #[case] expected: &str) {
    assert_eq!(expected, float_to_power_of_ten_with_error(input, error, digits));
  }
}
