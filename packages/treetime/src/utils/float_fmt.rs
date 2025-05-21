use lazy_static::lazy_static;
use pretty_dtoa::{FmtFloatConfig, dtoa};

lazy_static! {
  static ref FLOAT_CONFIG: FmtFloatConfig = FmtFloatConfig::default()
    .force_no_e_notation()
    .add_point_zero(true)
    .max_significant_digits(3)
    .radix_point('.')
    .round();
}

#[allow(clippy::string_slice)]
fn trim_trailing_zeros(input: String) -> String {
  match input.find('.') {
    Some(pos) => format!(
      "{}{}",
      &input[..pos],
      &input[pos..].trim_end_matches('0').trim_end_matches('.')
    ),
    None => input,
  }
}

fn float_format<F: Into<f64>>(x: F, config: FmtFloatConfig) -> String {
  trim_trailing_zeros(dtoa(x.into(), config))
}

pub fn float_to_significant_digits<F: Into<f64>>(x: F, max_significant_digits: u8) -> String {
  float_format(x, FLOAT_CONFIG.max_significant_digits(max_significant_digits))
}

pub fn float_to_decimal_digits<F: Into<f64>>(x: F, max_decimal_digits: i8) -> String {
  float_format(x, FLOAT_CONFIG.max_decimal_digits(max_decimal_digits))
}

pub fn float_to_digits<F: Into<f64>>(
  x: F,
  max_significant_digits: Option<u8>,
  max_decimal_digits: Option<i8>,
) -> String {
  let mut config = *FLOAT_CONFIG;
  if let Some(max_significant_digits) = max_significant_digits {
    config = config.max_significant_digits(max_significant_digits);
  }
  if let Some(max_decimal_digits) = max_decimal_digits {
    config = config.max_decimal_digits(max_decimal_digits);
  }
  float_format(x, config)
}
