#[cfg(test)]
mod tests {
  use crate::json_float::JsonFloat;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::finite(      1.5,               json!(1.5))]
  #[case::negative(    -2e-3,             json!(-2e-3))]
  #[case::pos_infinity(f64::INFINITY,     json!("inf"))]
  #[case::neg_infinity(f64::NEG_INFINITY, json!("-inf"))]
  #[trace]
  fn test_json_float_serializes_numbers_and_named_infinities(#[case] value: f64, #[case] expected: Value) {
    assert_eq!(expected, serde_json::to_value(JsonFloat(value)).unwrap());
  }

  #[test]
  fn test_json_float_serializes_nan_as_string() {
    assert_eq!(json!("nan"), serde_json::to_value(JsonFloat(f64::NAN)).unwrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::finite(      1.5,               json!(1.5))]
  #[case::pos_infinity(f64::INFINITY,     json!("inf"))]
  #[case::neg_infinity(f64::NEG_INFINITY, json!("-inf"))]
  #[trace]
  fn test_json_float_deserializes_what_it_serializes(#[case] expected: f64, #[case] json: Value) {
    assert_eq!(JsonFloat(expected), serde_json::from_value::<JsonFloat>(json).unwrap());
  }

  #[test]
  fn test_json_float_deserializes_nan() {
    assert!(serde_json::from_value::<JsonFloat>(json!("nan")).unwrap().0.is_nan());
  }

  #[test]
  fn test_json_float_rejects_other_strings() {
    assert_error!(
      serde_json::from_value::<JsonFloat>(json!("infinity")).map_err(Report::from),
      "expected a number, \"inf\", \"-inf\" or \"nan\", found \"infinity\""
    );
  }
}
