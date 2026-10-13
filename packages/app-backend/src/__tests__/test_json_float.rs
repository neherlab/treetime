#[cfg(test)]
mod tests {
  use crate::json_float::JsonFloat;
  use treetime_utils::io::json::from_json_value;
  use treetime_utils::io::json::to_json_value;

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
    assert_eq!(expected, to_json_value(&JsonFloat(value)).unwrap());
  }

  #[test]
  fn test_json_float_serializes_nan_as_string() {
    assert_eq!(json!("nan"), to_json_value(&JsonFloat(f64::NAN)).unwrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::finite(      1.5,               json!(1.5))]
  #[case::pos_infinity(f64::INFINITY,     json!("inf"))]
  #[case::neg_infinity(f64::NEG_INFINITY, json!("-inf"))]
  #[trace]
  fn test_json_float_deserializes_what_it_serializes(#[case] expected: f64, #[case] json: Value) {
    assert_eq!(JsonFloat(expected), from_json_value::<JsonFloat>(&json).unwrap());
  }

  #[test]
  fn test_json_float_deserializes_nan() {
    assert!(from_json_value::<JsonFloat>(&json!("nan")).unwrap().0.is_nan());
  }

  #[test]
  fn test_json_float_rejects_other_strings() {
    assert_error!(
      from_json_value::<JsonFloat>(&json!("infinity")),
      "When converting a JSON value: Unexpected: invalid value: expected a number, \"inf\", \"-inf\" or \"nan\", found \"infinity\""
    );
  }
}
