#[cfg(test)]
mod tests {
  use crate::runs::store::default_title;
  use chrono::{FixedOffset, TimeZone};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::utc(         0,         (2026, 9, 28, 14, 3, 59), "Run 2026-09-28 14:03")]
  #[case::east_of_utc( 2 * 3600,  (2026, 1,  2,  0, 0,  0), "Run 2026-01-02 00:00")]
  #[case::west_of_utc(-5 * 3600,  (2025, 12, 31, 23, 59, 0), "Run 2025-12-31 23:59")]
  #[trace]
  fn test_store_default_title_is_the_local_creation_minute(
    #[case] offset_seconds: i32,
    #[case] (year, month, day, hour, minute, second): (i32, u32, u32, u32, u32, u32),
    #[case] expected: &str,
  ) {
    let zone = FixedOffset::east_opt(offset_seconds).unwrap();
    let created_at = zone.with_ymd_and_hms(year, month, day, hour, minute, second).unwrap();
    assert_eq!(expected, default_title(&created_at));
  }
}
