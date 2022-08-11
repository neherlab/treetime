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

pub fn float_to_significant_digits<F: Into<f64>>(weight: F, max_significant_digits: u8) -> String {
  dtoa(
    weight.into(),
    FLOAT_CONFIG.max_significant_digits(max_significant_digits),
  )
}
