use num_traits::Float;
use pretty_dtoa::{FmtFloatConfig, dtoa};
use std::sync::LazyLock;

static FLOAT_CONFIG: LazyLock<FmtFloatConfig> =
  LazyLock::new(|| FmtFloatConfig::default().add_point_zero(true).radix_point('.').round());

pub trait FloatFormatExt {
  fn to_significant_digits(self, max_significant_digits: u8) -> String;

  fn to_decimal_digits(self, max_decimal_digits: i8) -> String;

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

pub fn float_to_digits<F: Into<f64>>(
  x: F,
  max_significant_digits: Option<u8>,
  max_decimal_digits: Option<i8>,
) -> String {
  let mut config = *FLOAT_CONFIG;

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

pub fn float_to_exponential<F: Into<f64>>(x: F, decimals: usize) -> String {
  let x = x.into();
  let formatted = format!("{x:.decimals$e}");
  let Some((mantissa, exponent)) = formatted.split_once('e') else {
    return formatted;
  };
  let (sign, digits) = exponent
    .strip_prefix('-')
    .map_or(('+', exponent), |digits| ('-', digits));
  format!("{mantissa}e{sign}{digits:0>2}")
}

pub fn float_to_power_of_ten(x: f64, significant_digits: usize) -> String {
  match power_of_ten_parts(x, significant_digits) {
    PowerOfTen::Scaled { mantissa, exponent } => with_power_of_ten(&signed(&mantissa), exponent),
    PowerOfTen::Plain(text) => text,
  }
}

pub fn float_to_power_of_ten_with_error(x: f64, error: f64, significant_digits: usize) -> String {
  match power_of_ten_parts(x, significant_digits) {
    PowerOfTen::Scaled { mantissa, exponent } => {
      let decimals = significant_digits.saturating_sub(1);
      let scaled_error = error / 10_f64.powi(exponent);
      with_power_of_ten(&format!("{} ± {scaled_error:.decimals$}", signed(&mantissa)), exponent)
    },
    PowerOfTen::Plain(text) => format!("{text} ± {}", float_to_power_of_ten(error, significant_digits)),
  }
}

enum PowerOfTen {
  Scaled { mantissa: String, exponent: i32 },
  Plain(String),
}

fn power_of_ten_parts(x: f64, significant_digits: usize) -> PowerOfTen {
  if x.is_nan() {
    return PowerOfTen::Plain("NaN".to_owned());
  }
  if x.is_infinite() {
    return PowerOfTen::Plain(signed(if x > 0.0 { "∞" } else { "-∞" }));
  }
  if x == 0.0 {
    return PowerOfTen::Plain("0".to_owned());
  }
  let decimals = significant_digits.saturating_sub(1);
  let formatted = format!("{x:.decimals$e}");
  match formatted
    .split_once('e')
    .map(|(mantissa, exponent)| (mantissa, exponent.parse::<i32>()))
  {
    Some((mantissa, Ok(exponent))) => PowerOfTen::Scaled {
      mantissa: mantissa.to_owned(),
      exponent,
    },
    _ => PowerOfTen::Plain(signed(&formatted)),
  }
}

fn with_power_of_ten(mantissa: &str, exponent: i32) -> String {
  if exponent == 0 {
    return mantissa.to_owned();
  }
  let superscript: String = exponent.to_string().chars().map(superscript_char).collect();
  format!("{mantissa} × 10{superscript}")
}

fn signed(text: &str) -> String {
  text.replace('-', "\u{2212}")
}

fn superscript_char(c: char) -> char {
  match c {
    '0' => '⁰',
    '1' => '¹',
    '2' => '²',
    '3' => '³',
    '4' => '⁴',
    '5' => '⁵',
    '6' => '⁶',
    '7' => '⁷',
    '8' => '⁸',
    '9' => '⁹',
    '-' => '⁻',
    other => other,
  }
}

fn float_format<F: Into<f64>>(x: F, config: FmtFloatConfig) -> String {
  let raw = dtoa(x.into(), config);
  trim_trailing_zeros(&raw)
}

#[expect(
  clippy::string_slice,
  reason = "the index comes from str::find on the ASCII exponent marker"
)]
pub(super) fn trim_trailing_zeros(input: &str) -> String {
  let (mantissa, exponent) = match input.find(['e', 'E']) {
    Some(pos) => (&input[..pos], Some(&input[pos..])),
    None => (input, None),
  };

  if exponent.is_some() {
    return input.to_owned();
  }

  match mantissa.find('.') {
    Some(pos) => {
      format!(
        "{}{}",
        &mantissa[..pos],
        mantissa[pos..].trim_end_matches('0').trim_end_matches('.')
      )
    },
    None => mantissa.to_owned(),
  }
}
