#[cfg(test)]
mod tests {
  use crate::config::catalog::{ListItemKind, RecordFieldSpec, SettingKind, SettingOption, value_spec};
  use app_commands::json_value::JsonValue;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_setting_records_fields_follow_the_item_schema() {
    let schema = json!({
      "properties": {
        "entries": {
          "type": "array",
          "items": { "$ref": "#/$defs/Entry" },
        },
      },
      "$defs": {
        "Entry": {
          "type": "object",
          "properties": {
            "name": { "type": "string", "description": "Name of the entry" },
            "sites": {
              "type": "array",
              "items": { "type": "string" },
              "examples": [["1-300"]],
              "x-unset": true,
            },
            "model": { "$ref": "#/$defs/Model", "x-unset": true, "default": "infer" },
            "weight": { "type": "integer", "description": "Weight of the entry\n\nA larger weight counts more." },
          },
          "required": ["name", "weight"],
          "additionalProperties": false,
        },
        "Model": {
          "oneOf": [
            { "const": "infer", "description": "Infer the model" },
            { "const": "jc69", "description": "Jukes-Cantor" },
          ],
        },
      },
    });

    let value = value_spec(&schema, &schema["properties"]["entries"], "entries").unwrap();

    assert_eq!(
      (
        SettingKind::Records,
        vec![
          RecordFieldSpec {
            key: o!("name"),
            label: o!("Name"),
            kind: SettingKind::Text,
            nullable: false,
            options: vec![],
            item_kind: ListItemKind::String,
            documented_default: None,
            examples: vec![],
            help: o!("Name of the entry"),
            more: o!(""),
          },
          RecordFieldSpec {
            key: o!("sites"),
            label: o!("Sites"),
            kind: SettingKind::List,
            nullable: true,
            options: vec![],
            item_kind: ListItemKind::String,
            documented_default: None,
            examples: vec![JsonValue(json!(["1-300"]))],
            help: o!(""),
            more: o!(""),
          },
          RecordFieldSpec {
            key: o!("model"),
            label: o!("Substitution model"),
            kind: SettingKind::Enum,
            nullable: true,
            options: vec![
              SettingOption {
                value: o!("infer"),
                help: o!("Infer the model"),
              },
              SettingOption {
                value: o!("jc69"),
                help: o!("Jukes-Cantor"),
              },
            ],
            item_kind: ListItemKind::String,
            documented_default: Some(JsonValue(json!("infer"))),
            examples: vec![],
            help: o!(""),
            more: o!(""),
          },
          RecordFieldSpec {
            key: o!("weight"),
            label: o!("Weight"),
            kind: SettingKind::Integer,
            nullable: false,
            options: vec![],
            item_kind: ListItemKind::String,
            documented_default: None,
            examples: vec![],
            help: o!("Weight of the entry"),
            more: o!("A larger weight counts more."),
          },
        ]
      ),
      (value.kind, value.fields)
    );
  }

  #[test]
  fn test_setting_records_field_not_in_required_can_be_left_out() {
    let schema = json!({
      "items": {
        "type": "object",
        "properties": { "flag": { "type": "boolean" } },
      },
      "type": "array",
    });

    let fields = value_spec(&schema, &schema, "entries").unwrap().fields;

    assert_eq!(
      vec![(SettingKind::Switch, true)],
      fields
        .iter()
        .map(|field| (field.kind, field.nullable))
        .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_setting_records_other_kinds_have_no_fields() {
    let schema = json!({ "type": "array", "items": { "type": "number" } });

    let value = value_spec(&schema, &schema, "relax").unwrap();

    assert_eq!(
      (SettingKind::List, ListItemKind::Number, vec![]),
      (value.kind, value.item_kind, value.fields)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::object_setting(
    json!({ "type": "object", "properties": { "name": { "type": "string" } } }),
    "setting `entries` has a schema construct the settings form does not render"
  )]
  #[case::items_without_fields(
    json!({ "type": "array", "items": { "type": "object", "additionalProperties": { "type": "string" } } }),
    "list setting `entries` has objects without named fields as items, which the settings form does not render"
  )]
  #[case::object_field(
    json!({ "type": "array", "items": { "type": "object", "properties": { "inner": { "type": "object", "properties": {} } } } }),
    "setting `entries.inner` has a schema construct the settings form does not render"
  )]
  #[case::records_field(
    json!({ "type": "array", "items": { "type": "object", "properties": {
      "inner": { "type": "array", "items": { "type": "object", "properties": { "name": { "type": "string" } } } },
    } } }),
    "field `entries.inner` is a list of objects inside a list of objects, which the settings form does not render"
  )]
  #[case::list_of_lists_field(
    json!({ "type": "array", "items": { "type": "object", "properties": {
      "inner": { "type": "array", "items": { "type": "array", "items": { "type": "string" } } },
    } } }),
    "list setting `entries.inner` has items the settings form does not render"
  )]
  #[trace]
  fn test_setting_records_rejects_shapes_the_form_does_not_render(#[case] schema: Value, #[case] message: &str) {
    assert_error!(value_spec(&schema, &schema, "entries"), message);
  }
}
