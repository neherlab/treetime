use lazy_static::lazy_static;
use maplit::btreemap;
use regex::{Error, Regex};
use std::collections::BTreeMap;

/// Try to convert strftime-style format string to a corresponding regular expression.
///
/// Note that at the current state the implementation is quite silly and not all formats are supported.
///
/// See: https://man7.org/linux/man-pages/man3/strftime.3.html
///
/// See: https://docs.rs/chrono/latest/chrono/format/strftime/index.html
pub fn date_format_to_regex(format: &str) -> Result<Regex, Error> {
  let escaped_format = regex::escape(format);
  let result = SPEC_MAP.iter().fold(escaped_format, |acc, (spec, expr)| {
    acc.replace(&format!("%{spec}"), expr)
  });
  Regex::new(&format!("^{result}$"))
}

lazy_static! {
  static ref SPEC_MAP: BTreeMap<char, &'static str> = btreemap! {
    'Y' => r"(?P<year>[\dX]{4})",
    'y' => r"(?P<year_short>[\dX]{2})",
    'm' => r"(?P<month>[\dX]{2})",
    'd' => r"(?P<day>[\dX]{2})",
    'H' => r"(?P<hour>[\dX]{2})",
    'M' => r"(?P<minute>[\dX]{2})",
    'S' => r"(?P<second>[\dX]{2})",
    'w' => r"(?P<week_day>[\dX]{1})",
    'W' => r"(?P<week>[\dX]{2})",
  };
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case("2024-07-23", "%Y-%m-%d", r"^(?P<year>[\dX]{4})\-(?P<month>[\dX]{2})\-(?P<day>[\dX]{2})$")]
  #[case("2024/07/23", "%Y/%m/%d", r"^(?P<year>[\dX]{4})/(?P<month>[\dX]{2})/(?P<day>[\dX]{2})$"  )]
  #[case("2024.07.23", "%Y.%m.%d", r"^(?P<year>[\dX]{4})\.(?P<month>[\dX]{2})\.(?P<day>[\dX]{2})$")]
  #[case("23-07-2024", "%d-%m-%Y", r"^(?P<day>[\dX]{2})\-(?P<month>[\dX]{2})\-(?P<year>[\dX]{4})$")]
  #[case("23/07/2024", "%d/%m/%Y", r"^(?P<day>[\dX]{2})/(?P<month>[\dX]{2})/(?P<year>[\dX]{4})$"  )]
  #[case("23.07.2024", "%d.%m.%Y", r"^(?P<day>[\dX]{2})\.(?P<month>[\dX]{2})\.(?P<year>[\dX]{4})$")]
  #[case("20240723"  , "%Y%m%d",   r"^(?P<year>[\dX]{4})(?P<month>[\dX]{2})(?P<day>[\dX]{2})$"    )]
  #[trace]
  fn test_date_format_to_regex(#[case] example: &str, #[case] fmt: &str, #[case] expected: &str) {
    let regex = date_format_to_regex(fmt).unwrap();
    assert_eq!(expected, regex.as_str());

    assert!(regex.is_match(example));

    let caps = regex.captures(example).unwrap();
    let year = caps.name("year").map_or("None", |s| s.as_str());
    let month = caps.name("month").map_or("None", |s| s.as_str());
    let day = caps.name("day").map_or("None", |s| s.as_str());
    let actual= format!("{year}-{month}-{day}");
    assert_eq!("2024-07-23", actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case("2024-07-XX", "%Y-%m-%d", r"^(?P<year>[\dX]{4})\-(?P<month>[\dX]{2})\-(?P<day>[\dX]{2})$")]
  #[case("2024-XX-XX", "%Y-%m-%d", r"^(?P<year>[\dX]{4})\-(?P<month>[\dX]{2})\-(?P<day>[\dX]{2})$")]
  #[trace]
  fn test_date_format_to_regex_uncertain(#[case] example: &str, #[case] fmt: &str, #[case] expected: &str) {
    let regex = date_format_to_regex(fmt).unwrap();
    assert_eq!(expected, regex.as_str());

    assert!(regex.is_match(example));

    let caps = regex.captures(example).unwrap();
    let year = caps.name("year").map_or("None", |s| s.as_str());
    let month = caps.name("month").map_or("None", |s| s.as_str());
    let day = caps.name("day").map_or("None", |s| s.as_str());
    let actual= format!("{year}-{month}-{day}");
    assert_eq!(example, actual);
  }
}
