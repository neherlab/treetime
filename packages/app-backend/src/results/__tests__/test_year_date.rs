#[cfg(test)]
mod tests {
  use crate::results::year_date::YearDate;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::start_of_year(     2016.0,            "2016-01-01")]
  #[case::half_a_leap_year(  2016.5,            "2016-07-02")]
  #[case::noon_of_a_day(     2015.0 + 171.5 / 365.0, "2015-06-21")]
  #[trace]
  fn test_year_date_calendar_day(#[case] year: f64, #[case] expected: &str) {
    assert_eq!(YearDate { year, date: o!(expected) }, YearDate::new(year));
  }
}
