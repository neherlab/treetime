#[cfg(test)]
mod tests {
  use crate::datetime::datetime::iso;
  use crate::datetime::year_fraction::{
    date_to_year_fraction, year_fraction_at_month_start, year_fraction_days_between, year_fraction_to_date,
    year_fraction_to_datestring,
  };
  use crate::{assert_error, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};
  use chrono::{Datelike, NaiveDate, SubsecRound};
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::regression_year_before_zero_day_150_of_365(     -1909.59,                  "-01910-05-30" )]
  #[case::year_zero_is_a_leap_year_day_184_of_366(           0.5,                    "0000-07-02"   )]
  #[case::whole_negative_year(                              -5.0,                    "-00005-01-01" )]
  #[case::just_below_zero_is_the_last_day_of_year_minus_one( -0.001,                 "-00001-12-31" )]
  #[case::leap_day_of_year_minus_4(                         -4.0   + 59.5 / 366.0,   "-00004-02-29" )]
  #[case::leap_day_of_year_minus_400(                     -400.0   + 59.5 / 366.0,   "-00400-02-29" )]
  #[case::no_leap_day_in_year_minus_100(                  -100.0   + 59.5 / 365.0,   "-00100-03-01" )]
  #[case::deep_past_day_256_of_365(                      -5869.0   + 255.5 / 365.0,  "-05869-09-13" )]
  #[case::last_day_of_four_digit_years(                   9999.999,                  "9999-12-31"   )]
  #[case::first_five_digit_year_is_a_leap_year(          10000.0   + 183.5 / 366.0,  "+10000-07-02" )]
  #[case::first_year_of_chrono(                        -262143.0   + 182.5 / 365.0,  "-262143-07-02")]
  #[case::last_year_of_chrono(                          262142.0   + 182.5 / 365.0,  "+262142-07-02")]
  #[case::fraction_rounds_up_to_the_next_year(               -1e-20,                 "0000-01-01"   )]
  #[case::four_digit_year(                                2012.0   + 166.5 / 366.0,  "2012-06-15"   )]
  #[trace]
  fn test_year_fraction_to_datestring(#[case] year_fraction: f64, #[case] expected: &str) {
    assert_eq!(Some(expected.to_owned()), year_fraction_to_datestring(year_fraction));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nan(                      f64::NAN         )]
  #[case::infinity(                 f64::INFINITY    )]
  #[case::negative_infinity(        f64::NEG_INFINITY)]
  #[case::year_before_chrono(  -262143.5            )]
  #[case::year_after_chrono(    262143.0            )]
  #[case::huge(                    1e300            )]
  #[trace]
  fn test_year_fraction_to_datestring_without_calendar_date(#[case] year_fraction: f64) {
    assert_eq!(None, year_fraction_to_datestring(year_fraction));
  }

  #[test]
  fn test_year_fraction_to_date_covers_the_range_of_chrono() {
    let first = date_to_year_fraction(&noon(NaiveDate::MIN));
    let last = date_to_year_fraction(&noon(NaiveDate::MAX));
    assert_eq!(
      (Some(noon(NaiveDate::MIN)), Some(noon(NaiveDate::MAX))),
      (
        year_fraction_to_date(first).map(|date| date.round_subsecs(0)),
        year_fraction_to_date(last).map(|date| date.round_subsecs(0))
      )
    );
    assert_eq!(
      (None, None),
      (
        year_fraction_to_date(f64::from(NaiveDate::MIN.year()) - 1e-6),
        year_fraction_to_date(f64::from(NaiveDate::MAX.year()) + 1.0)
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::first_day_of_year_one(       (     1,  1,  1),      1.0 +   0.5 / 365.0)]
  #[case::day_before_skipped_leap_day( (  1900,  2, 28),   1900.0 +  58.5 / 365.0)]
  #[case::day_after_skipped_leap_day(  (  1900,  3,  1),   1900.0 +  59.5 / 365.0)]
  #[case::last_day_of_common_year(     (  1900, 12, 31),   1900.0 + 364.5 / 365.0)]
  #[case::leap_day_of_400_year_cycle(  (  2000,  2, 29),   2000.0 +  59.5 / 366.0)]
  #[case::day_after_leap_day(          (  2000,  3,  1),   2000.0 +  60.5 / 366.0)]
  #[case::last_day_of_leap_year(       (  2000, 12, 31),   2000.0 + 365.5 / 366.0)]
  #[case::middle_of_common_year(       (  2023,  7,  2),   2023.0 + 182.5 / 365.0)]
  #[case::july_2_of_leap_year(         (  2024,  7,  2),   2024.0 + 183.5 / 366.0)]
  #[case::first_day_of_year_zero(      (     0,  1,  1),      0.0 +   0.5 / 366.0)]
  #[case::last_day_of_year_zero(       (     0, 12, 31),      0.0 + 365.5 / 366.0)]
  #[case::last_day_of_year_minus_one(  (    -1, 12, 31),     -1.0 + 364.5 / 365.0)]
  #[case::leap_day_of_year_minus_4(    (    -4,  2, 29),     -4.0 +  59.5 / 366.0)]
  #[case::day_60_of_year_minus_100(    (  -100,  3,  1),   -100.0 +  59.5 / 365.0)]
  #[case::first_day_of_chrono(         (NaiveDate::MIN.year(),  1,  1), f64::from(NaiveDate::MIN.year()) +   0.5 / 365.0)]
  #[case::last_day_of_chrono(          (NaiveDate::MAX.year(), 12, 31), f64::from(NaiveDate::MAX.year()) + 364.5 / 365.0)]
  #[trace]
  fn test_date_to_year_fraction_is_noon_of_the_day(#[case] (year, month, day): (i32, u32, u32), #[case] expected: f64) {
    pretty_assert_ulps_eq!(expected, date_to_year_fraction(&noon(ymd(year, month, day))), max_ulps = 4);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::first_day_of_year_one(       1.0 +   0.5 / 365.0, (     1,  1,  1))]
  #[case::day_after_skipped_leap_day(  1900.0 +  59.5 / 365.0, (  1900,  3,  1))]
  #[case::leap_day_of_400_year_cycle(  2000.0 +  59.5 / 366.0, (  2000,  2, 29))]
  #[case::last_day_of_leap_year(       2000.0 + 365.5 / 366.0, (  2000, 12, 31))]
  #[case::first_day_of_year_zero(         0.0 +   0.5 / 366.0, (     0,  1,  1))]
  #[case::last_day_of_year_minus_one(    -1.0 + 364.5 / 365.0, (    -1, 12, 31))]
  #[case::leap_day_of_year_minus_4(      -4.0 +  59.5 / 366.0, (    -4,  2, 29))]
  #[trace]
  fn test_year_fraction_to_date_at_noon(#[case] year_fraction: f64, #[case] (year, month, day): (i32, u32, u32)) {
    assert_eq!(
      Some(noon(ymd(year, month, day))),
      year_fraction_to_date(year_fraction).map(|date| date.round_subsecs(0))
    );
  }

  proptest! {
    #[test]
    fn test_prop_year_fraction_noon_roundtrip(
      day in NaiveDate::MIN.num_days_from_ce()..=NaiveDate::MAX.num_days_from_ce()
    ) {
      let date = noon(NaiveDate::from_num_days_from_ce_opt(day).unwrap());
      let actual = year_fraction_to_date(date_to_year_fraction(&date)).map(|date| date.round_subsecs(0));
      prop_assert_eq!(Some(date), actual);
    }
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::across_year_zero(       (-1.0 + 364.5 / 365.0,   0.0 + 0.5 / 366.0),           1.0)]
  #[case::across_leap_day(        (-4.0 +  58.5 / 366.0,  -4.0 + 60.5 / 366.0),          2.0)]
  #[case::across_leap_day_back(   (-4.0 +  60.5 / 366.0,  -4.0 + 58.5 / 366.0),         -2.0)]
  #[case::thousand_cycles_of_400( (-200_000.0,             200_000.0),       146_097_000.0)]
  #[case::thousand_cycles_back(   ( 200_000.0,            -200_000.0),      -146_097_000.0)]
  #[trace]
  fn test_year_fraction_days_between(#[case] (from, to): (f64, f64), #[case] expected: f64) {
    pretty_assert_abs_diff_eq!(expected, year_fraction_days_between(from, to).unwrap(), epsilon = 1e-13);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::to_without_date(  (0.0,      1e300))]
  #[case::from_without_date((f64::NAN, 0.0  ))]
  #[trace]
  fn test_year_fraction_days_between_without_date(#[case] (from, to): (f64, f64)) {
    assert_eq!(None, year_fraction_days_between(from, to));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::january_is_the_year(     (2016,  1), 2016.0)]
  #[case::march_after_a_leap_day(  (2024,  3), 2024.0 + 60.0 / 366.0)]
  #[case::march_without_a_leap_day((2023,  3), 2023.0 + 59.0 / 365.0)]
  #[case::december(                (2023, 12), 2023.0 + 334.0 / 365.0)]
  #[case::year_before_the_era(     (-4,    7), -4.0 + 182.0 / 366.0)]
  #[trace]
  fn test_year_fraction_at_month_start(#[case] (year, month): (i32, u32), #[case] expected: f64) {
    pretty_assert_ulps_eq!(expected, year_fraction_at_month_start(year, month).unwrap(), max_ulps = 4);
  }

  #[test]
  fn test_year_fraction_at_month_start_is_the_first_instant_of_the_month() {
    assert_eq!(
      Some(iso("2016-09-01T00:00:00.000Z")),
      year_fraction_to_date(year_fraction_at_month_start(2016, 9).unwrap()).map(|date| date.round_subsecs(0))
    );
  }

  #[test]
  fn test_year_fraction_at_month_start_rejects_month_thirteen() {
    assert_error!(
      year_fraction_at_month_start(2016, 13),
      "Month 13 of year 2016 is not a calendar month"
    );
  }

  mod helpers {
    use chrono::{DateTime, NaiveDate, TimeZone, Utc};

    pub(super) fn ymd(year: i32, month: u32, day: u32) -> NaiveDate {
      NaiveDate::from_ymd_opt(year, month, day).unwrap()
    }

    pub(super) fn noon(date: NaiveDate) -> DateTime<Utc> {
      Utc.from_utc_datetime(&date.and_hms_opt(12, 0, 0).unwrap())
    }
  }

  use helpers::{noon, ymd};
}
