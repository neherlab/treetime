use crate::datetime::date_range::DateRange;
use crate::datetime::options::DateParserOptions;
use crate::datetime::parse_uncertain_date::parse_date_uncertain;
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
