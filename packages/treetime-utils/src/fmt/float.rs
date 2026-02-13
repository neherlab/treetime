use lazy_static::lazy_static;
use num_traits::Float;
use pretty_dtoa::{FmtFloatConfig, dtoa};

lazy_static! {
  static ref FLOAT_CONFIG: FmtFloatConfig = FmtFloatConfig::default().add_point_zero(true).radix_point('.').round();
}

#[allow(clippy::string_slice)]
fn trim_trailing_zeros(input: &str) -> String {
  let (mantissa, exponent) = match input.find(['e', 'E']) {
    Some(pos) => (&input[..pos], Some(&input[pos..])),
    None => (input, None),
  };

  // Never trim trailing zeros in scientific notation.
  if exponent.is_some() {
    return input.to_owned();
  }

  match mantissa.find('.') {
    Some(pos) => {
      // For plain decimals we can remove the fractional part entirely.
      format!(
        "{}{}",
        &mantissa[..pos],
        &mantissa[pos..].trim_end_matches('0').trim_end_matches('.')
      )
    },
    None => mantissa.to_owned(),
  }
}

fn float_format<F: Into<f64>>(x: F, config: FmtFloatConfig) -> String {
  let raw = dtoa(x.into(), config);
  trim_trailing_zeros(&raw)
}

/// Trait providing convenient float formatting methods
pub trait FloatFormatExt {
  /// Format float to a specific number of significant digits
  fn to_significant_digits(self, max_significant_digits: u8) -> String;

  /// Format float to a specific number of decimal digits
  fn to_decimal_digits(self, max_decimal_digits: i8) -> String;

  /// Format float with optional significant and decimal digit limits
  fn to_digits(self, max_significant_digits: Option<u8>, max_decimal_digits: Option<i8>) -> String;
}

impl<F: Float + Into<f64>> FloatFormatExt for F {
  fn to_significant_digits(self, max_significant_digits: u8) -> String {
    float_to_digits(self, Some(max_significant_digits), None)
  }

  fn to_decimal_digits(self, max_decimal_digits: i8) -> String {
    float_to_digits(self, None, Some(max_decimal_digits))
  }

  fn to_digits(self, max_significant_digits: Option<u8>, max_decimal_digits: Option<i8>) -> String {
    float_to_digits(self, max_significant_digits, max_decimal_digits)
  }
}

pub fn float_to_significant_digits<F: Into<f64>>(x: F, max_significant_digits: u8) -> String {
  float_to_digits(x, Some(max_significant_digits), None)
}

pub fn float_to_decimal_digits<F: Into<f64>>(x: F, max_decimal_digits: i8) -> String {
  float_to_digits(x, None, Some(max_decimal_digits))
}

pub fn float_to_digits<F: Into<f64>>(
  x: F,
  max_significant_digits: Option<u8>,
  max_decimal_digits: Option<i8>,
) -> String {
  let mut config = *FLOAT_CONFIG;

  // If neither constraint is specified, use the default significant digits from FLOAT_CONFIG
  if max_significant_digits.is_none() && max_decimal_digits.is_none() {
    config = config.max_significant_digits(3);
  }

  if let Some(max_significant_digits) = max_significant_digits {
    config = config.max_significant_digits(max_significant_digits);
  }

  if let Some(max_decimal_digits) = max_decimal_digits {
    config = config.max_decimal_digits(max_decimal_digits);
  }

  float_format(x, config)
}

#[cfg(test)]
mod tests {
  use super::*;
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
