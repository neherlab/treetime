use crate::datetime::date_range::DateRange;
use crate::datetime::datetime::days_in_month;
use crate::datetime::format_to_regex::date_format_to_regex;
use crate::datetime::options::DateParserOptions;
use crate::make_error;
use eyre::Report;
use itertools::Itertools;
use lazy_static::lazy_static;
use regex::Regex;

/// Try to read date with uncertain components, e.g. 2024-07-XX
pub fn parse_date_uncertain(date_uncertain_str: &str, _options: &DateParserOptions) -> Result<DateRange, Report> {
  for regex in DATE_UNCERTAIN_REGEXES.iter() {
    if let Some(caps) = regex.captures(date_uncertain_str) {
      let (year, month, day) = ["year", "month", "day"]
        .iter()
        .map(|name| {
          caps.name(name).and_then(|s| {
            let s = s.as_str();
            if s.chars().all(|c| c == 'X') { None } else { Some(s) }
          })
        })
        .collect_tuple()
        .unwrap();

      if let Ok(Some(range)) = determine_date_range(year, month, day) {
        return Ok(range);
      }
    }
  }
  make_error!("Unrecognized date format: {date_uncertain_str}")
}

fn determine_date_range(
  year: Option<&str>,
  month: Option<&str>,
  day: Option<&str>,
) -> Result<Option<DateRange>, Report> {
  let year = year.and_then(|year| resolve_uncertain_date_component(year, (1, 9999)).ok());

  let month = month.and_then(|month| resolve_uncertain_date_component(month, (1, 12)).ok());

  let day = if let (Some(year), Some(month)) = (year, month) {
    let day_max = days_in_month(year.1, month.1)?;
    day.and_then(|day| resolve_uncertain_date_component(day, (1, day_max)).ok())
  } else {
    None
  };

  match (year, month, day) {
    (Some((year_min, year_max)), Some((month_min, month_max)), Some((day_min, day_max))) => Ok(Some(
      DateRange::from_ymd((year_min, month_min, day_min), (year_max, month_max, day_max)),
    )),
    (Some((year_min, year_max)), Some((month_min, month_max)), None) => {
      let day_min = 1;
      let day_max = days_in_month(year_max, month_max).ok();
      if let Some(day_max) = day_max {
        Ok(Some(DateRange::from_ymd(
          (year_min, month_min, day_min),
          (year_max, month_max, day_max),
        )))
      } else {
        Ok(None)
      }
    },
    (Some((year_min, year_max)), None, None) => Ok(Some(DateRange::from_ymd((year_min, 1, 1), (year_max, 12, 31)))),
    _ => Ok(None),
  }
}

/// Try to resolve uncertainty of the given date component string and clamp it to the given bounds.
/// For example, month `XX` should be clamped to (1, 12).
fn resolve_uncertain_date_component(s: impl AsRef<str>, bounds: (u32, u32)) -> Result<(u32, u32), Report> {
  let s = s.as_ref();
  let min = s.replace('X', "0").parse::<u32>()?.clamp(bounds.0, bounds.1);
  let max = s.replace('X', "9").parse::<u32>()?.clamp(bounds.0, bounds.1);
  Ok((min, max))
}

lazy_static! {
  static ref DATE_UNCERTAIN_REGEXES: Vec<Regex> = create_date_uncertain_regexes();
}

fn create_date_uncertain_regexes() -> Vec<Regex> {
  #[rustfmt::skip]
  const FORMATS: &[&str] = &[
    "%Y-%m-%d",
    "%Y/%m/%d",
    "%Y.%m.%d",
    //
    "%Y-%m",
    "%Y/%m",
    "%Y.%m",
    //
    "%d-%m-%Y",
    "%d/%m/%Y",
    "%d.%m.%Y",
    //
    "%m-%Y",
    "%m/%Y",
    "%m.%Y",
    //
    "%Y",
    //
    "%Y%m%d",
  ];

  FORMATS.iter().map(|f| date_format_to_regex(f).unwrap()).collect_vec()
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::datetime::date_range::DateRange;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  fn r(begin: &str, end: &str) -> DateRange {
    DateRange::from_iso(begin, end)
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::ymd_dash_full("2024-07-23", r("2024-07-23T00:00:00.000Z", "2024-07-23T23:59:59.999999999Z"))]
  #[case::ymd_dash_day_xx("2024-07-XX", r("2024-07-01T00:00:00.000Z", "2024-07-31T23:59:59.999999999Z"))]
  #[case::ymd_dash_month_day_xx("2024-XX-XX", r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  #[case::ymd_dash_decade_xx("20XX-XX-XX", r("2000-01-01T00:00:00.000Z", "2099-12-31T23:59:59.999999999Z"))]
  #[case::ymd_dash_millennium_xx("2XXX-XX-XX", r("2000-01-01T00:00:00.000Z", "2999-12-31T23:59:59.999999999Z"))]
  //
  #[case::ym_dash_full("2024-07",    r("2024-07-01T00:00:00.000Z", "2024-07-31T23:59:59.999999999Z"))]
  #[case::ym_dash_month_xx("2024-XX",    r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  //
  #[case::ymd_slash_full("2024/07/23", r("2024-07-23T00:00:00.000Z", "2024-07-23T23:59:59.999999999Z"))]
  #[case::ymd_slash_day_xx("2024/07/XX", r("2024-07-01T00:00:00.000Z", "2024-07-31T23:59:59.999999999Z"))]
  #[case::ymd_slash_month_day_xx("2024/XX/XX", r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  #[case::ymd_slash_decade_xx("20XX/XX/XX", r("2000-01-01T00:00:00.000Z", "2099-12-31T23:59:59.999999999Z"))]
  #[case::ymd_slash_millennium_xx("2XXX/XX/XX", r("2000-01-01T00:00:00.000Z", "2999-12-31T23:59:59.999999999Z"))]
  //
  #[case::ym_slash_full("2024/07",    r("2024-07-01T00:00:00.000Z", "2024-07-31T23:59:59.999999999Z"))]
  #[case::ym_slash_month_xx("2024/XX",    r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  //
  #[case::ymd_dot_full("2024.07.23", r("2024-07-23T00:00:00.000Z", "2024-07-23T23:59:59.999999999Z"))]
  #[case::ymd_dot_day_xx("2024.07.XX", r("2024-07-01T00:00:00.000Z", "2024-07-31T23:59:59.999999999Z"))]
  #[case::ymd_dot_month_day_xx("2024.XX.XX", r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  #[case::ymd_dot_decade_xx("20XX.XX.XX", r("2000-01-01T00:00:00.000Z", "2099-12-31T23:59:59.999999999Z"))]
  #[case::ymd_dot_millennium_xx("2XXX.XX.XX", r("2000-01-01T00:00:00.000Z", "2999-12-31T23:59:59.999999999Z"))]
  //
  #[case::ym_dot_full("2024.07",    r("2024-07-01T00:00:00.000Z", "2024-07-31T23:59:59.999999999Z"))]
  #[case::ym_dot_month_xx("2024.XX",    r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  //
  #[case::dmy_dash_full("23-07-2024", r("2024-07-23T00:00:00.000Z", "2024-07-23T23:59:59.999999999Z"))]
  #[case::dmy_dash_day_xx("XX-07-2024", r("2024-07-01T00:00:00.000Z", "2024-07-31T23:59:59.999999999Z"))]
  #[case::dmy_dash_day_month_xx("XX-XX-2024", r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  #[case::dmy_dash_decade_xx("XX-XX-20XX", r("2000-01-01T00:00:00.000Z", "2099-12-31T23:59:59.999999999Z"))]
  #[case::dmy_dash_millennium_xx("XX-XX-2XXX", r("2000-01-01T00:00:00.000Z", "2999-12-31T23:59:59.999999999Z"))]
  //
  #[case::my_dash_full("07-2024",    r("2024-07-01T00:00:00.000Z", "2024-07-31T23:59:59.999999999Z"))]
  #[case::my_dash_month_xx("XX-2024",    r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  //
  #[case::dmy_slash_full("23/07/2024", r("2024-07-23T00:00:00.000Z", "2024-07-23T23:59:59.999999999Z"))]
  #[case::dmy_slash_day_xx("XX/07/2024", r("2024-07-01T00:00:00.000Z", "2024-07-31T23:59:59.999999999Z"))]
  #[case::dmy_slash_day_month_xx("XX/XX/2024", r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  #[case::dmy_slash_decade_xx("XX/XX/20XX", r("2000-01-01T00:00:00.000Z", "2099-12-31T23:59:59.999999999Z"))]
  #[case::dmy_slash_millennium_xx("XX/XX/2XXX", r("2000-01-01T00:00:00.000Z", "2999-12-31T23:59:59.999999999Z"))]
  //
  #[case::my_slash_full("07/2024",    r("2024-07-01T00:00:00.000Z", "2024-07-31T23:59:59.999999999Z"))]
  #[case::my_slash_month_xx("XX/2024",    r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  //
  #[case::dmy_dot_full("23.07.2024", r("2024-07-23T00:00:00.000Z", "2024-07-23T23:59:59.999999999Z"))]
  #[case::dmy_dot_day_xx("XX.07.2024", r("2024-07-01T00:00:00.000Z", "2024-07-31T23:59:59.999999999Z"))]
  #[case::dmy_dot_day_month_xx("XX.XX.2024", r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  #[case::dmy_dot_decade_xx("XX.XX.20XX", r("2000-01-01T00:00:00.000Z", "2099-12-31T23:59:59.999999999Z"))]
  #[case::dmy_dot_millennium_xx("XX.XX.2XXX", r("2000-01-01T00:00:00.000Z", "2999-12-31T23:59:59.999999999Z"))]
  //
  #[case::my_dot_full("07.2024",    r("2024-07-01T00:00:00.000Z", "2024-07-31T23:59:59.999999999Z"))]
  #[case::my_dot_month_xx("XX.2024",    r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  //
  #[case::year_only("2024",       r("2024-01-01T00:00:00.000Z", "2024-12-31T23:59:59.999999999Z"))]
  #[trace]
  fn test_date_parse_uncertain(#[case] input: &str, #[case] expected: DateRange) {
    let options = DateParserOptions::default();
    let actual = parse_date_uncertain(input, &options).unwrap();
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::all_x_ymd_dash("XXXX-XX-XX")]
  #[case::all_x_ymd_slash("XXXX/XX/XX")]
  #[case::all_x_ymd_dot("XXXX.XX.XX")]
  #[case::all_x_ym_dash("XXXX-XX")]
  #[case::all_x_ym_slash("XXXX/XX")]
  #[case::all_x_ym_dot("XXXX.XX")]
  #[case::all_x_dmy_dash("XX-XX-XXXX")]
  #[case::all_x_dmy_slash("XX/XX/XXXX")]
  #[case::all_x_dmy_dot("XX.XX.XXXX")]
  #[case::all_x_my_dash("XX-XXXX")]
  #[case::all_x_my_slash("XX/XXXX")]
  #[case::all_x_my_dot("XX.XXXX")]
  #[case::all_x_year("XXXX")]
  #[case::empty("")]
  #[case::whitespace("  ")]
  #[trace]
  fn test_date_parse_uncertain_error(#[case] input: &str) {
    let options = DateParserOptions::default();
    let actual = parse_date_uncertain(input, &options).ok();
    assert_eq!(None, actual);
  }
}
