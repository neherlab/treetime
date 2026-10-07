use crate::error::{NewickWriteError, WriteContext, write_error};
use pretty_dtoa::{FmtFloatConfig, dtoa};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NumberFormat {
  pub significant_digits: Option<u8>,
  pub decimal_digits: Option<i8>,
  pub point_zero: bool,
}

impl NumberFormat {
  pub fn format(&self, value: f64) -> Result<String, NewickWriteError> {
    format_number(value, *self)
  }
}

pub(crate) fn format_shortest(value: f64) -> Result<String, NewickWriteError> {
  format_number(value, NumberFormat::default())
}

fn format_number(value: f64, format: NumberFormat) -> Result<String, NewickWriteError> {
  if !value.is_finite() {
    return Err(write_error!("Newick cannot represent the number {value}"));
  }
  let mut config = FmtFloatConfig::default()
    .add_point_zero(format.point_zero)
    .radix_point('.');
  if let Some(significant_digits) = format.significant_digits {
    if significant_digits == 0 {
      return Err(write_error!("The number of significant digits must be at least 1"));
    }
    config = config.max_significant_digits(significant_digits);
  }
  let value = match format.decimal_digits {
    Some(decimal_digits) => round_to_decimal_digits(value, decimal_digits)?,
    None => value,
  };
  let text = dtoa(value, config);
  Ok(trim_fraction_zeros(text, format.point_zero))
}

fn trim_fraction_zeros(text: String, point_zero: bool) -> String {
  if text.contains(['e', 'E']) || !text.contains('.') {
    return text;
  }
  let trimmed = text.trim_end_matches('0');
  match (trimmed.strip_suffix('.'), point_zero) {
    (Some(_), true) => format!("{trimmed}0"),
    (Some(integer), false) => integer.to_owned(),
    (None, _) => trimmed.to_owned(),
  }
}

fn round_to_decimal_digits(value: f64, decimal_digits: i8) -> Result<f64, NewickWriteError> {
  if decimal_digits >= 0 {
    let decimals = usize::from(decimal_digits.unsigned_abs());
    let rounded = format!("{value:.decimals$}");
    return rounded
      .parse::<f64>()
      .with_context(|| format!("When rounding {value} to {decimals} decimal digits"));
  }
  let scale = 10_f64.powi(i32::from(decimal_digits.unsigned_abs()));
  Ok((value / scale).round() * scale)
}
