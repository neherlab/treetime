use crate::utils::datetime::date_range::DateRange;
use chrono::{DateTime, Datelike, TimeZone, Utc};
use chronoutil::RelativeDuration;
use std::time::Duration as StdDuration;
use time::util::days_in_year;

/// Convert DateTime object to a year-fraction number
///
/// NOTE: the calculation is not reciprocal to `date_to_year_fraction()` due to precision loss in
/// floating point arithmetics. In practical applications the precision is ensured only up to seconds.
pub fn date_to_year_fraction(date: &DateTime<Utc>) -> f64 {
  let year = date.year();
  assert!(year >= 1); // TODO: implement BC dates?
  let frac = (date.ordinal() as f64 - 0.5) / days_in_year(year) as f64;
  year as f64 + frac
}

pub fn date_range_to_year_fraction_range(date_range: &DateRange) -> (f64, f64) {
  let begin = date_to_year_fraction(date_range.begin());
  let end = date_to_year_fraction(date_range.end());
  (begin, end)
}

/// Convert year-fraction to DateTime object
///
/// NOTE: the calculation is not reciprocal to `date_to_year_fraction()` due to precision loss in
/// floating point arithmetics. In practical applications the precision is ensured only up to seconds.
pub fn year_fraction_to_date(year_fraction: f64) -> DateTime<Utc> {
  let year = year_fraction.trunc() as i32;
  let fraction = year_fraction.fract();
  let seconds_in_year = (days_in_year(year) as u64) * 24 * 60 * 60;
  let seconds_since_year_start = seconds_in_year as f64 * fraction;
  let dt = RelativeDuration::from(StdDuration::from_secs_f64(seconds_since_year_start));
  Utc.ymd(year, 1, 1).and_hms(0, 0, 0) + dt
}

#[cfg(test)]
mod tests {
  #![allow(clippy::excessive_precision)]
  use super::*;
  use crate::pretty_assert_ulps_eq;
  use crate::utils::datetime::datetime::iso;
  use chrono::{SubsecRound, Utc};
  use pretty_assertions::assert_eq;
  use rstest::*;

  #[rustfmt::skip]
  #[rstest]
  #[case::start_of_ad                                 (iso("0001-01-01T12:00:00.000Z"),     1.0013698630)]
  #[case::pompeii_destruction                         (iso("0079-08-24T12:00:00.000Z"),    79.6452054794)]
  #[case::founding_of_baghdad                         (iso("0762-07-30T12:00:00.000Z"),   762.5767123287)]
  #[case::supernova_1054_ad                           (iso("1054-07-04T12:00:00.000Z"),  1054.5054794520)]
  #[case::first_crusade                               (iso("1096-08-15T12:00:00.000Z"),  1096.6215846994)]
  #[case::gregorian_calendar_introduced               (iso("1582-10-15T12:00:00.000Z"),  1582.7876712328)]
  #[case::completion_of_the_suez_canal                (iso("1869-11-17T12:00:00.000Z"),  1869.8780821917)]
  #[case::wright_brothers_first_flight                (iso("1903-12-17T12:00:00.000Z"),  1903.9602739726)]
  #[case::discovery_of_penicillin                     (iso("1928-09-28T12:00:00.000Z"),  1928.7418032786)]
  #[case::invention_of_the_transistor                 (iso("1947-12-23T12:00:00.000Z"),  1947.9767123287)]
  #[case::intel_microprocessor_release                (iso("1971-09-06T12:00:00.000Z"),  1971.6808219178)]
  #[case::first_iss_launch                            (iso("1998-11-20T12:00:00.000Z"),  1998.8863013698)]
  #[case::dot_com_bubble_burst                        (iso("2000-03-10T12:00:00.000Z"),  2000.1898907103)]
  #[case::iphone_launch                               (iso("2007-06-29T12:00:00.000Z"),  2007.4917808219)]
  #[case::this_test_written                           (iso("2024-09-05T12:00:00.000Z"),  2024.6789617486)]
  #[case::halley_comet_return                         (iso("2061-07-28T12:00:00.000Z"),  2061.5712328767)]
  #[case::transit_of_venus                            (iso("2117-12-11T12:00:00.000Z"),  2117.9438356164)]
  #[case::pluto_full_orbit_since_discovery            (iso("2178-03-23T12:00:00.000Z"),  2178.2232876712)]
  #[case::thousand_years_since_this_test              (iso("3024-09-05T12:00:00.000Z"),  3024.6789617486)]
  //
  #[case::half_century                                (iso("1950-07-02T12:00:00.000Z"),  1950.5000000000)]
  #[case::quarter_century                             (iso("2025-04-01T12:00:00.000Z"),  2025.2479452054)]
  #[case::three_quarters_century                      (iso("2075-09-30T12:00:00.000Z"),  2075.7465753424)]
  #[case::hapl_millennium                             (iso("2500-12-31T12:00:00.000Z"),  2500.9986301369)]
  #[case::halfway_non_leap_year                       (iso("2023-07-02T12:00:00.000Z"),  2023.5000000000)]
  //
  #[case::leap_year_2024_day_before                   (iso("2024-02-28T12:00:00.000Z"),  2024.15983606553)]
  #[case::leap_year_2024_leap_day                     (iso("2024-02-29T12:00:00.000Z"),  2024.16256830600)]
  #[case::leap_year_2024_day_after                    (iso("2024-03-01T12:00:00.000Z"),  2024.16530054648)]
  #[case::leap_year_2024_halfway                      (iso("2024-07-02T12:00:00.000Z"),  2024.50136612028)]
  #[case::leap_year_2024_last_day                     (iso("2024-12-31T12:00:00.000Z"),  2024.99863387971)]
  #[case::leap_year_2024_day_after_year               (iso("2025-01-01T12:00:00.000Z"),  2025.00136986303)]
  //
  #[case::round_leap_year_2000_day_before             (iso("2000-02-28T12:00:00.000Z"),  2000.15983606557)]
  #[case::round_leap_year_2000_leap_day               (iso("2000-02-29T12:00:00.000Z"),  2000.16256830601)]
  #[case::round_leap_year_2000_day_after              (iso("2000-03-01T12:00:00.000Z"),  2000.16530054644)]
  #[case::round_leap_year_2000_halfway                (iso("2000-07-02T12:00:00.000Z"),  2000.50136612021)]
  #[case::round_leap_year_2000_last_day               (iso("2000-12-31T12:00:00.000Z"),  2000.99863387978)]
  #[case::round_leap_year_2000_day_after_year         (iso("2001-01-01T12:00:00.000Z"),  2001.00136986301)]
  //
  #[case::excluded_leap_year_1900_day_before          (iso("1900-02-27T12:00:00.000Z"),  1900.15753424657)]
  #[case::excluded_leap_year_1900_want_to_be_leap_day (iso("1900-02-28T12:00:00.000Z"),  1900.16027397260)]
  #[case::excluded_leap_year_1900_day_after           (iso("1900-03-01T12:00:00.000Z"),  1900.16301369863)]
  #[case::excluded_leap_year_1900_last_day            (iso("1900-12-31T12:00:00.000Z"),  1900.99863013698)]
  #[case::excluded_leap_year_1900_halfway             (iso("1900-07-02T12:00:00.000Z"),  1900.50000000000)]
  #[case::excluded_leap_year_1900_day_after_year      (iso("1901-01-01T12:00:00.000Z"),  1901.00136986301)]
  //
  #[case::non_leap_year_2023_day_before               (iso("2023-02-27T12:00:00.000Z"),  2023.15753424657)]
  #[case::non_leap_year_2023_want_to_be_leap_day      (iso("2023-02-28T12:00:00.000Z"),  2023.16027397260)]
  #[case::non_leap_year_2023_day_after                (iso("2023-03-01T12:00:00.000Z"),  2023.16301369863)]
  #[case::non_leap_year_2023_halfway                  (iso("2023-07-02T12:00:00.000Z"),  2023.50000000000)]
  #[case::non_leap_year_2023_last_day                 (iso("2023-12-31T12:00:00.000Z"),  2023.99863013698)]
  #[case::non_leap_year_2023_day_after_year           (iso("2024-01-01T12:00:00.000Z"),  2024.00136612021)]
  //
  #[case::y2k_day_before                              (iso("1999-12-31T12:00:00.000Z"),  1999.99863013698)]
  #[case::y2k                                         (iso("2000-01-01T12:00:00.000Z"),  2000.00136612021)]
  #[case::y2k_day_after                               (iso("2000-01-02T12:00:00.000Z"),  2000.00409836065)]
  #[case::y2k38_before                                (iso("2038-01-17T12:00:00.000Z"),  2038.04520547945)]
  #[case::y2k38                                       (iso("2038-01-18T12:00:00.000Z"),  2038.04794520547)]
  #[case::y2k38_after                                 (iso("2038-01-19T12:00:00.000Z"),  2038.05068493150)]
  //
  #[case::unix_epoch_day_before                       (iso("1969-12-31T12:00:00.000Z"),  1969.99863013698)]
  #[case::unix_epoch                                  (iso("1970-01-01T12:00:00.000Z"),  1970.00136986301)]
  #[case::unix_epoch_day_after                        (iso("1970-01-02T12:00:00.000Z"),  1970.00410958904)]
  //
  #[trace]
  fn test_date_to_year_fraction(#[case] date: DateTime<Utc>, #[case] yf: f64) {
    let actual = date_to_year_fraction(&date);
    pretty_assert_ulps_eq!(yf, actual, epsilon = 1e-6);

    let actual = year_fraction_to_date(yf).round_subsecs(0); // Rounded to seconds; precision loss
    assert_eq!(date, actual);
  }
}
