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
