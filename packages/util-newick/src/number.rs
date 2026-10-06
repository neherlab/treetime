use eyre::{Report, WrapErr, eyre};
use pretty_dtoa::{FmtFloatConfig, dtoa};

pub(crate) fn format_number(
  value: f64,
  significant_digits: Option<u8>,
  decimal_digits: Option<i8>,
) -> Result<String, Report> {
  if !value.is_finite() {
    return Err(eyre!("Newick cannot represent the number {value}"));
  }
  let mut config = FmtFloatConfig::default().add_point_zero(false).radix_point('.');
  if let Some(significant_digits) = significant_digits {
    if significant_digits == 0 {
      return Err(eyre!("The number of significant digits must be at least 1"));
    }
    config = config.max_significant_digits(significant_digits);
  }
  let value = match decimal_digits {
    Some(decimal_digits) => round_to_decimal_digits(value, decimal_digits)?,
    None => value,
  };
  Ok(dtoa(value, config))
}

fn round_to_decimal_digits(value: f64, decimal_digits: i8) -> Result<f64, Report> {
  if decimal_digits >= 0 {
    let decimals = usize::from(decimal_digits.unsigned_abs());
    let rounded = format!("{value:.decimals$}");
    return rounded
      .parse::<f64>()
      .wrap_err_with(|| format!("When rounding {value} to {decimals} decimal digits"));
  }
  let scale = 10_f64.powi(i32::from(decimal_digits.unsigned_abs()));
  Ok((value / scale).round() * scale)
}

pub(crate) fn format_shortest(value: f64) -> Result<String, Report> {
  format_number(value, None, None)
}

#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(
    result_defaulted,
    reason = "text that is not a finite number is a string value, not a malformed number"
  )
)]
pub(crate) fn parse_number_text(text: &str) -> Option<f64> {
  let is_number_syntax = text
    .chars()
    .all(|c| c.is_ascii_digit() || matches!(c, '+' | '-' | '.' | 'e' | 'E'))
    && text.chars().any(|c| c.is_ascii_digit());
  if !is_number_syntax {
    return None;
  }
  text.parse::<f64>().ok().filter(|value| value.is_finite())
}
