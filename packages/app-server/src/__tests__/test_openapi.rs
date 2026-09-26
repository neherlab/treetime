#[cfg(test)]
mod tests {
  use crate::openapi::add_discriminators;
  use crate::routes::api_doc;
  use helpers::{
    discriminated, keyword_locations, response_components, schemas_after, tagged_unions_without_discriminator,
  };
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use treetime_utils::assert_error;

  #[test]
  fn test_openapi_every_tagged_union_has_a_discriminator() {
    let doc = api_doc().unwrap();
    assert_eq!(Vec::<String>::new(), tagged_unions_without_discriminator(&doc));
  }

  #[test]
  fn test_openapi_discriminated_unions_are_the_tagged_serde_enums() {
    let doc = api_doc().unwrap();
    assert_eq!(
      vec![
        ("CheckConfigResponse", "status"),
        ("CoalescentPrior", "kind"),
        ("CommandResults", "command"),
        ("JobEvent", "type"),
        ("OperationRequest", "operation"),
        ("RunConfigResponse", "status"),
        ("RunEvent", "type"),
        ("SettingDifference", "kind"),
        ("TerminalEvent", "status"),
      ],
      discriminated(&doc)
    );
  }

  #[test]
  fn test_openapi_run_event_variants_carry_the_fields_of_the_event() {
    let doc = api_doc().unwrap();
    let schemas = &doc["components"]["schemas"];
    assert_eq!(
      (
        json!({
          "propertyName": "type",
          "mapping": {
            "started": "#/components/schemas/RunEventStarted",
            "progress": "#/components/schemas/RunEventProgress",
            "log": "#/components/schemas/RunEventLog",
            "iteration": "#/components/schemas/RunEventIteration",
            "terminal": "#/components/schemas/RunEventTerminal",
          },
        }),
        json!(["seq", "time", "type", "data"]),
        json!(["seq", "time", "type", "data"]),
      ),
      (
        schemas["RunEvent"]["discriminator"].clone(),
        schemas["RunEventTerminal"]["required"].clone(),
        Value::Array(
          schemas["RunEventTerminal"]["properties"]
            .as_object()
            .unwrap()
            .keys()
            .map(|key| json!(key))
            .collect()
        ),
      )
    );
  }

  #[test]
  fn test_openapi_response_schemas_carry_no_defaults() {
    let doc = api_doc().unwrap();
    let components = response_components(&doc);
    let with_defaults = components
      .iter()
      .flat_map(|name| keyword_locations(&doc["components"]["schemas"][name.as_str()], "default", name))
      .collect::<Vec<_>>();
    assert_eq!(
      (true, true, Vec::<String>::new()),
      (
        components.contains("RunRecord"),
        components.contains("RunEventTerminal"),
        with_defaults
      )
    );
  }

  #[test]
  fn test_openapi_request_schemas_carry_defaults() {
    let doc = api_doc().unwrap();
    assert_eq!(
      vec![
        "CreateRunRequest/properties/title/default",
        "CreateRunRequest/properties/defer_start/default"
      ],
      keyword_locations(
        &doc["components"]["schemas"]["CreateRunRequest"],
        "default",
        "CreateRunRequest"
      )
    );
  }

  #[test]
  fn test_openapi_internally_tagged_union_becomes_named_variants() {
    let actual = schemas_after(json!({
      "Shape": {
        "description": "A shape.",
        "oneOf": [
          {
            "description": "A circle.",
            "type": "object",
            "properties": { "r": { "type": "number" }, "kind": { "type": "string", "const": "circle" } },
            "required": ["kind", "r"],
          },
          {
            "type": "object",
            "properties": { "kind": { "type": "string", "const": "rounded-square" } },
            "required": ["kind"],
            "additionalProperties": false,
          },
        ],
      },
    }))
    .unwrap();
    assert_eq!(
      json!({
        "Shape": {
          "description": "A shape.",
          "oneOf": [
            { "$ref": "#/components/schemas/ShapeCircle" },
            { "$ref": "#/components/schemas/ShapeRoundedSquare" },
          ],
          "discriminator": {
            "propertyName": "kind",
            "mapping": {
              "circle": "#/components/schemas/ShapeCircle",
              "rounded-square": "#/components/schemas/ShapeRoundedSquare",
            },
          },
        },
        "ShapeCircle": {
          "description": "A circle.",
          "type": "object",
          "properties": { "r": { "type": "number" }, "kind": { "type": "string", "const": "circle" } },
          "required": ["kind", "r"],
        },
        "ShapeRoundedSquare": {
          "type": "object",
          "properties": { "kind": { "type": "string", "const": "rounded-square" } },
          "required": ["kind"],
          "additionalProperties": false,
        },
      }),
      actual
    );
  }

  #[test]
  fn test_openapi_fields_next_to_a_flattened_union_move_into_each_variant() {
    let actual = schemas_after(json!({
      "Event": {
        "description": "An event.",
        "type": "object",
        "properties": { "seq": { "type": "integer" } },
        "required": ["seq"],
        "oneOf": [
          {
            "type": "object",
            "properties": { "type": { "type": "string", "const": "a" }, "seq": { "type": "integer" } },
            "required": ["type", "seq"],
          },
          {
            "type": "object",
            "properties": { "type": { "type": "string", "const": "b" }, "data": { "type": "string" } },
            "required": ["type", "data"],
          },
        ],
      },
    }))
    .unwrap();
    assert_eq!(
      (
        json!({
          "type": "object",
          "properties": { "seq": { "type": "integer" }, "type": { "type": "string", "const": "a" } },
          "required": ["seq", "type"],
        }),
        json!({
          "type": "object",
          "properties": {
            "seq": { "type": "integer" },
            "type": { "type": "string", "const": "b" },
            "data": { "type": "string" },
          },
          "required": ["seq", "type", "data"],
        }),
        json!({
          "description": "An event.",
          "oneOf": [{ "$ref": "#/components/schemas/EventA" }, { "$ref": "#/components/schemas/EventB" }],
          "discriminator": {
            "propertyName": "type",
            "mapping": { "a": "#/components/schemas/EventA", "b": "#/components/schemas/EventB" },
          },
        }),
      ),
      (
        actual["EventA"].clone(),
        actual["EventB"].clone(),
        actual["Event"].clone()
      )
    );
  }

  #[test]
  fn test_openapi_unions_without_one_tag_stay_unchanged() {
    let schemas = json!({
      "Level": {
        "oneOf": [{ "type": "string", "const": "low" }, { "type": "string", "const": "high" }],
      },
      "Untagged": {
        "oneOf": [
          { "type": "object", "properties": { "a": { "type": "string" } }, "required": ["a"] },
          { "type": "object", "properties": { "b": { "type": "string" } }, "required": ["b"] },
        ],
      },
      "OptionalTag": {
        "oneOf": [
          { "type": "object", "properties": { "kind": { "const": "a" } }, "required": ["kind"] },
          { "type": "object", "properties": { "kind": { "const": "b" } } },
        ],
      },
      "Holder": {
        "type": "object",
        "properties": { "level": { "anyOf": [{ "$ref": "#/components/schemas/Level" }, { "type": "null" }] } },
      },
    });
    assert_eq!(schemas.clone(), schemas_after(schemas).unwrap());
  }

  #[test]
  fn test_openapi_union_with_two_tag_properties_is_an_error() {
    assert_error!(
      schemas_after(json!({
        "Pair": {
          "oneOf": [
            {
              "type": "object",
              "properties": { "a": { "const": "x" }, "b": { "const": "y" } },
              "required": ["a", "b"],
            },
            {
              "type": "object",
              "properties": { "a": { "const": "z" }, "b": { "const": "w" } },
              "required": ["a", "b"],
            },
          ],
        },
      })),
      "the union `Pair` can be told apart by each of `a`, `b`; the discriminator needs exactly one tag property"
    );
  }

  #[test]
  fn test_openapi_union_with_a_repeated_tag_value_is_an_error() {
    assert_error!(
      schemas_after(json!({
        "Twice": {
          "oneOf": [
            { "type": "object", "properties": { "kind": { "const": "a" } }, "required": ["kind"] },
            {
              "type": "object",
              "properties": { "kind": { "const": "a" }, "n": { "type": "integer" } },
              "required": ["kind", "n"],
            },
          ],
        },
      })),
      "the tagged union `Twice` has two variants with the same `kind`"
    );
  }

  #[test]
  fn test_openapi_variant_named_like_another_schema_is_an_error() {
    assert_error!(
      schemas_after(json!({
        "Shape": {
          "oneOf": [
            { "type": "object", "properties": { "kind": { "const": "circle" } }, "required": ["kind"] },
            { "type": "object", "properties": { "kind": { "const": "square" } }, "required": ["kind"] },
          ],
        },
        "ShapeCircle": { "type": "string" },
      })),
      "the variant `ShapeCircle` of the tagged union `Shape` has the name of another schema; rename one of the types"
    );
  }

  #[test]
  fn test_openapi_union_with_an_unknown_keyword_next_to_one_of_is_an_error() {
    assert_error!(
      schemas_after(json!({
        "Open": {
          "additionalProperties": false,
          "oneOf": [
            { "type": "object", "properties": { "kind": { "const": "a" } }, "required": ["kind"] },
            { "type": "object", "properties": { "kind": { "const": "b" } }, "required": ["kind"] },
          ],
        },
      })),
      "the tagged union `Open` has the keyword `additionalProperties` next to `oneOf`, which its variants cannot take over"
    );
  }

  #[test]
  fn test_openapi_shared_field_described_differently_by_a_variant_is_an_error() {
    assert_error!(
      schemas_after(json!({
        "Event": {
          "type": "object",
          "properties": { "seq": { "type": "integer" } },
          "required": ["seq"],
          "oneOf": [
            {
              "type": "object",
              "properties": { "type": { "const": "a" }, "seq": { "type": "string" } },
              "required": ["type"],
            },
            { "type": "object", "properties": { "type": { "const": "b" } }, "required": ["type"] },
          ],
        },
      })),
      "the tagged union `Event` describes the field `seq` twice, differently"
    );
  }

  #[test]
  fn test_openapi_shared_fields_next_to_a_closed_variant_is_an_error() {
    assert_error!(
      schemas_after(json!({
        "Event": {
          "type": "object",
          "properties": { "seq": { "type": "integer" } },
          "required": ["seq"],
          "oneOf": [
            {
              "type": "object",
              "properties": { "type": { "const": "a" } },
              "required": ["type"],
              "additionalProperties": false,
            },
            { "type": "object", "properties": { "type": { "const": "b" } }, "required": ["type"] },
          ],
        },
      })),
      "the tagged union `Event` has fields next to a variant that allows no other fields"
    );
  }

  mod helpers {
    use super::add_discriminators;
    use aide::openapi::{Components, OpenApi};
    use eyre::Report;
    use itertools::Itertools;
    use serde_json::{Map, Value, json};
    use std::collections::BTreeSet;

    const COMPONENTS_PREFIX: &str = "#/components/schemas/";

    pub(super) fn schemas_after(schemas: Value) -> Result<Value, Report> {
      let mut api = OpenApi {
        components: Some(Components {
          schemas: serde_json::from_value(schemas)?,
          ..Components::default()
        }),
        ..OpenApi::default()
      };
      add_discriminators(&mut api)?;
      Ok(serde_json::to_value(api)?["components"]["schemas"].clone())
    }

    pub(super) fn discriminated(doc: &Value) -> Vec<(&str, &str)> {
      doc["components"]["schemas"]
        .as_object()
        .unwrap()
        .iter()
        .filter_map(|(name, schema)| Some((name.as_str(), schema["discriminator"]["propertyName"].as_str()?)))
        .sorted()
        .collect()
    }

    pub(super) fn response_components(doc: &Value) -> BTreeSet<String> {
      let mut pending = doc["paths"]
        .as_object()
        .unwrap()
        .values()
        .flat_map(|item| item.as_object().unwrap().values())
        .flat_map(|operation| operation["responses"].as_object().into_iter().flatten())
        .flat_map(|(_, response)| references(response))
        .collect_vec();
      let mut seen = BTreeSet::new();
      while let Some(name) = pending.pop() {
        if seen.insert(name.clone()) {
          pending.extend(references(&doc["components"]["schemas"][name.as_str()]));
        }
      }
      seen
    }

    pub(super) fn keyword_locations(schema: &Value, keyword: &str, pointer: &str) -> Vec<String> {
      match schema {
        Value::Object(object) => object
          .iter()
          .flat_map(|(key, child)| {
            let location = format!("{pointer}/{key}");
            let here = (key == keyword && !pointer.ends_with("/properties")).then(|| location.clone());
            here.into_iter().chain(keyword_locations(child, keyword, &location))
          })
          .collect(),
        Value::Array(items) => items
          .iter()
          .enumerate()
          .flat_map(|(index, child)| keyword_locations(child, keyword, &format!("{pointer}/{index}")))
          .collect(),
        _ => vec![],
      }
    }

    fn references(value: &Value) -> Vec<String> {
      match value {
        Value::Object(object) => object
          .iter()
          .flat_map(|(key, child)| match (key.as_str(), child.as_str()) {
            ("$ref", Some(reference)) => reference
              .strip_prefix(COMPONENTS_PREFIX)
              .map(ToOwned::to_owned)
              .into_iter()
              .collect(),
            _ => references(child),
          })
          .collect(),
        Value::Array(items) => items.iter().flat_map(references).collect(),
        _ => vec![],
      }
    }

    pub(super) fn tagged_unions_without_discriminator(doc: &Value) -> Vec<String> {
      let schemas = doc["components"]["schemas"].as_object().unwrap();
      let mut found = vec![];
      collect_unions(doc, "#", schemas, &mut found);
      found
    }

    fn collect_unions(value: &Value, pointer: &str, schemas: &Map<String, Value>, found: &mut Vec<String>) {
      match value {
        Value::Object(object) => {
          if let Some(Value::Array(one_of)) = object.get("oneOf")
            && let Some(problem) = union_problem(object, one_of, schemas)
          {
            found.push(format!("{pointer}: {problem}"));
          }
          for (key, child) in object {
            collect_unions(child, &format!("{pointer}/{key}"), schemas, found);
          }
        },
        Value::Array(items) => {
          for (index, child) in items.iter().enumerate() {
            collect_unions(child, &format!("{pointer}/{index}"), schemas, found);
          }
        },
        _ => {},
      }
    }

    fn union_problem(union: &Map<String, Value>, one_of: &[Value], schemas: &Map<String, Value>) -> Option<String> {
      let variants = one_of
        .iter()
        .map(|variant| match variant["$ref"].as_str() {
          Some(reference) => schemas.get(reference.strip_prefix(COMPONENTS_PREFIX)?),
          None => Some(variant),
        })
        .collect::<Option<Vec<_>>>()?;
      let first = variants.first()?["properties"].as_object()?;
      let tag = first.keys().find(|property| {
        variants.iter().all(|variant| {
          variant["properties"][property.as_str()]["const"].is_string()
            && variant["required"]
              .as_array()
              .is_some_and(|required| required.contains(&json!(property)))
        })
      })?;
      let expected = variants
        .iter()
        .zip(one_of)
        .map(|(variant, item)| {
          (
            variant["properties"][tag.as_str()]["const"]
              .as_str()
              .unwrap()
              .to_owned(),
            item["$ref"].clone(),
          )
        })
        .collect_vec();
      let discriminator = union.get("discriminator").unwrap_or(&Value::Null);
      let mapping = discriminator["mapping"].as_object();
      let complete = discriminator["propertyName"] == json!(tag)
        && mapping.is_some_and(|mapping| {
          mapping.len() == expected.len()
            && expected
              .iter()
              .all(|(value, reference)| reference.is_string() && mapping.get(value) == Some(reference))
        });
      (!complete).then(|| format!("union tagged by `{tag}` has no complete discriminator"))
    }
  }
}
