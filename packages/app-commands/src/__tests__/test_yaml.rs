#[cfg(test)]
mod tests {
  use crate::yaml::yaml_document;
  use generators::json_value;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use serde_json::{Value, json};

  #[test]
  fn test_yaml_document_writes_nested_mappings_and_sequences() {
    let actual = yaml_document(&json!({
      "workspace": "/data/runs",
      "ui": { "theme": "dark", "sidebar_width": 400, "tags": ["a", "b"], "empty": {} },
    }))
    .unwrap();
    let expected = indoc! {r#"
      workspace: "/data/runs"
      ui:
        theme: "dark"
        sidebar_width: 400
        tags:
          - "a"
          - "b"
        empty: {}
    "#};
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::boolean_word(json!({ "true": "yes" }))]
  #[case::null_word(   json!({ "null": "~" }))]
  #[case::number_key(  json!({ "1": 1 }))]
  #[case::empty_key(   json!({ "": "" }))]
  #[case::colon_key(   json!({ "a: b": "c: d" }))]
  #[trace]
  fn test_yaml_document_reads_back_keys_that_look_like_other_scalars(#[case] value: Value) {
    let text = yaml_document(&value).unwrap();
    let actual: Value = serde_saphyr::from_str(&text).unwrap();
    assert_eq!(value, actual);
  }

  proptest! {
    #[test]
    fn test_prop_yaml_document_roundtrip(value in json_value()) {
      let text = yaml_document(&value).unwrap();
      let actual: Value = serde_saphyr::from_str(&text).unwrap();
      prop_assert_eq!(value, actual);
    }
  }

  mod generators {
    use proptest::prelude::*;
    use serde_json::{Map, Number, Value};

    pub(super) fn json_value() -> impl Strategy<Value = Value> {
      let leaf = prop_oneof![
        Just(Value::Null),
        any::<bool>().prop_map(Value::Bool),
        any::<i64>().prop_map(Value::from),
        any::<u64>().prop_map(Value::from),
        (-1e300_f64..1e300_f64).prop_filter_map("finite", |x| Number::from_f64(x).map(Value::Number)),
        any::<String>().prop_map(Value::String),
      ];
      leaf.prop_recursive(4, 32, 6, |inner| {
        prop_oneof![
          prop::collection::vec(inner.clone(), 0..6).prop_map(Value::Array),
          prop::collection::btree_map("[a-z_][a-z0-9_]{0,8}", inner, 0..6)
            .prop_map(|entries| Value::Object(entries.into_iter().collect::<Map<String, Value>>())),
        ]
      })
    }
  }
}
