use lazy_static::lazy_static;
use pretty_dtoa::{dtoa, FmtFloatConfig};

lazy_static! {
  static ref FLOAT_CONFIG: FmtFloatConfig = FmtFloatConfig::default()
    .force_no_e_notation()
    .add_point_zero(true)
    .max_significant_digits(3)
    .radix_point('.')
    .round();
}

fn float_format<F: Into<f64>>(x: F, config: FmtFloatConfig) -> String {
  dtoa(x.into(), config).trim_end_matches('0').to_owned()
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
