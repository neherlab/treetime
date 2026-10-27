#[cfg(test)]
mod tests {
  use crate::results::year_date::YearDate;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::start_of_year(          2016.0,                 Some("2016-01-01")  )]
  #[case::half_a_leap_year(       2016.5,                 Some("2016-07-02")  )]
  #[case::noon_of_a_day(          2015.0 + 171.5 / 365.0, Some("2015-06-21")  )]
  #[case::before_year_zero(      -1909.59,                Some("-01910-05-30"))]
  #[case::outside_calendar_range( 1e9,                    None                )]
  #[trace]
  fn test_year_date_calendar_day(#[case] year: f64, #[case] expected: Option<&str>) {
    assert_eq!(YearDate { year, date: expected.map(str::to_owned) }, YearDate::new(year));
  }
}
