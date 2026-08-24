#[cfg(test)]
mod __tests__;

use std::sync::LazyLock;

use num_traits::Float;
use pretty_dtoa::{FmtFloatConfig, dtoa};

static FLOAT_CONFIG: LazyLock<FmtFloatConfig> =
  LazyLock::new(|| FmtFloatConfig::default().add_point_zero(true).radix_point('.').round());

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
