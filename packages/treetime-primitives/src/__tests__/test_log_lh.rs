#[cfg(test)]
mod tests {
  use crate::LogLh;
  use std::mem::{align_of, size_of};
  use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_log_lh_construction_and_extraction() {
    let log_lh = LogLh::new(-12.5);

    pretty_assert_ulps_eq!(-12.5, log_lh.value(), max_ulps = 1);
  }

  #[test]
  fn test_log_lh_constants() {
    pretty_assert_ulps_eq!(0.0, LogLh::ZERO.value(), max_ulps = 1);
    assert!(LogLh::IMPOSSIBLE.value().is_infinite());
    assert!(LogLh::IMPOSSIBLE.value().is_sign_negative());
  }

  #[test]
  fn test_log_lh_addition_and_assignment() {
    let mut actual = LogLh::new(-2.0) + LogLh::new(-3.5);
    actual += LogLh::new(1.0);

    pretty_assert_ulps_eq!(-4.5, actual.value(), max_ulps = 1);
  }

  #[test]
  fn test_log_lh_owned_and_borrowed_sums() {
    let values = [LogLh::new(-1.0), LogLh::new(-2.0), LogLh::new(-3.0)];

    let owned: LogLh = values.into_iter().sum();
    let borrowed: LogLh = values.iter().sum();

    pretty_assert_ulps_eq!(-6.0, owned.value(), max_ulps = 1);
    pretty_assert_ulps_eq!(-6.0, borrowed.value(), max_ulps = 1);
  }

  #[test]
  fn test_log_lh_subtraction_returns_difference() {
    let difference = LogLh::new(-2.0) - LogLh::new(-5.5);

    pretty_assert_ulps_eq!(3.5, difference, max_ulps = 1);
  }

  #[test]
  fn test_log_lh_negation_returns_optimizer_cost() {
    let cost = -LogLh::new(-2.5);

    pretty_assert_ulps_eq!(2.5, cost, max_ulps = 1);
  }

  #[test]
  fn test_log_lh_comparison() {
    assert!(LogLh::new(-2.0) > LogLh::new(-3.0));
  }

  #[test]
  fn test_log_lh_layout_matches_f64() {
    assert_eq!(size_of::<f64>(), size_of::<LogLh>());
    assert_eq!(align_of::<f64>(), align_of::<LogLh>());
  }

  #[test]
  fn test_log_lh_json_is_transparent() {
    let expected = "-12.5";

    let actual = json_write_str(&LogLh::new(-12.5), JsonPretty(false)).expect("LogLh serialization should succeed");
    let roundtrip: LogLh = json_read_str(&actual).expect("LogLh deserialization should succeed");

    assert_eq!(expected, actual);
    pretty_assert_ulps_eq!(-12.5, roundtrip.value(), max_ulps = 1);
  }
}
