#[cfg(test)]
mod tests {
  use crate::routes::{api_doc, api_doc_with};
  use app_commands::config::schema::draft2020_settings;
  use helpers::{enum_spelling_problems, keyword_locations_in, open_schemas, schema_nodes, unreachable_components};
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};

  const OPEN_SCHEMAS_ALLOWED: &[&str] = &[
    "#/components/schemas/AuspiceDocument",
    "#/paths/~1api~1openapi.json/get/responses/200/content/application~1json/schema",
  ];

  #[test]
  fn test_openapi_contract_no_schema_is_open() {
    let doc = api_doc().unwrap();
    let open = open_schemas(&doc)
      .into_iter()
      .filter(|pointer| !OPEN_SCHEMAS_ALLOWED.iter().any(|allowed| pointer.starts_with(allowed)))
      .collect::<Vec<_>>();
    assert_eq!(Vec::<String>::new(), open);
  }

  #[test]
  fn test_openapi_contract_every_component_is_reachable_from_a_route_or_the_catalog() {
    let doc = api_doc().unwrap();
    assert_eq!(Vec::<String>::new(), unreachable_components(&doc, &["SettingCatalog"]));
  }

  #[test]
  fn test_openapi_contract_no_schema_allows_null() {
    let doc = api_doc().unwrap();
    let nullable = schema_nodes(&doc)
      .into_iter()
      .filter(|(_, node)| {
        let type_allows_null = match node.get("type") {
          Some(Value::String(single)) => single == "null",
          Some(Value::Array(types)) => types.contains(&json!("null")),
          _ => false,
        };
        let enum_allows_null = node
          .get("enum")
          .and_then(Value::as_array)
          .is_some_and(|values| values.contains(&Value::Null));
        type_allows_null
          || enum_allows_null
          || node.get("const") == Some(&Value::Null)
          || node.get("nullable") == Some(&json!(true))
      })
      .map(|(pointer, _)| pointer)
      .collect::<Vec<_>>();
    assert_eq!(Vec::<String>::new(), nullable);
  }

  #[test]
  fn test_openapi_contract_enums_share_one_spelling() {
    let doc = api_doc().unwrap();
    assert_eq!(Vec::<String>::new(), enum_spelling_problems(&doc));
  }

  #[test]
  fn test_openapi_contract_no_property_defaults_to_null() {
    let doc = api_doc().unwrap();
    assert_eq!(
      Vec::<String>::new(),
      keyword_locations_in(&doc, "default", &Value::Null)
    );
  }

  #[test]
  fn test_openapi_contract_no_optional_field_is_serialized_as_null() {
    let doc = api_doc_with(&draft2020_settings().for_serialize()).unwrap();
    let required_unset = schema_nodes(&doc)
      .into_iter()
      .flat_map(|(pointer, node)| {
        let required = node
          .get("required")
          .and_then(Value::as_array)
          .cloned()
          .unwrap_or_default();
        let properties = node
          .get("properties")
          .and_then(Value::as_object)
          .cloned()
          .unwrap_or_default();
        properties
          .into_iter()
          .filter(move |(name, property)| {
            required.contains(&json!(name)) && property.get("x-unset") == Some(&json!(true))
          })
          .map(move |(name, _)| format!("{pointer}/properties/{name}"))
      })
      .collect::<Vec<_>>();
    assert_eq!(Vec::<String>::new(), required_unset);
  }

  #[test]
  fn test_openapi_contract_open_schema_is_found() {
    let doc = json!({
      "components": { "schemas": {
        "Closed": { "type": "object", "properties": { "a": { "type": "string" } } },
        "Any": { "description": "Anything." },
        "Map": { "type": "object", "additionalProperties": true },
        "Nested": { "type": "object", "properties": { "value": {} } },
      } },
      "paths": {},
    });
    assert_eq!(
      vec![
        "#/components/schemas/Any",
        "#/components/schemas/Map",
        "#/components/schemas/Nested/properties/value",
      ],
      open_schemas(&doc)
    );
  }

  #[test]
  fn test_openapi_contract_enum_spelling_problems_are_found() {
    let doc = json!({
      "components": { "schemas": {
        "Kebab": { "type": "string", "enum": ["mat-pb", "nwk"] },
        "Pascal": { "type": "string", "enum": ["MatPb", "Nwk"] },
        "Mixed": { "oneOf": [{ "const": "MatPb" }, { "const": "nwk" }] },
        "Snake": { "type": "string", "enum": ["invalid_request", "not_found"] },
      } },
      "paths": {},
    });
    assert_eq!(
      vec![
        "#/components/schemas/Kebab and #/components/schemas/Mixed and #/components/schemas/Pascal spell the same \
         values differently",
        "#/components/schemas/Mixed mixes case styles",
      ],
      enum_spelling_problems(&doc)
    );
  }

  #[test]
  fn test_openapi_contract_unreachable_component_is_found() {
    let doc = json!({
      "components": { "schemas": {
        "Used": { "type": "object", "properties": { "child": { "$ref": "#/components/schemas/Child" } } },
        "Child": { "type": "string" },
        "Orphan": { "type": "string" },
      } },
      "paths": { "/api/x": { "get": { "responses": { "200": { "content": { "application/json": {
        "schema": { "$ref": "#/components/schemas/Used" }
      } } } } } } },
    });
    assert_eq!(vec!["Orphan"], unreachable_components(&doc, &[]));
  }

  mod helpers {
    use itertools::Itertools;
    use serde_json::Value;
    use std::collections::{BTreeMap, BTreeSet};

    const COMPONENTS_PREFIX: &str = "#/components/schemas/";

    const ANNOTATIONS: &[&str] = &[
      "description",
      "title",
      "examples",
      "default",
      "deprecated",
      "readOnly",
      "writeOnly",
    ];

    const SCHEMA_MAPS: &[&str] = &["properties", "patternProperties", "$defs", "dependentSchemas"];

    const SCHEMA_LISTS: &[&str] = &["anyOf", "oneOf", "allOf", "prefixItems"];

    const SCHEMA_SINGLES: &[&str] = &[
      "items",
      "additionalProperties",
      "not",
      "if",
      "then",
      "else",
      "contains",
      "propertyNames",
      "unevaluatedItems",
      "unevaluatedProperties",
    ];

    const METHODS: &[&str] = &["get", "put", "post", "delete", "patch", "head", "options"];

    pub(super) fn schema_nodes(doc: &Value) -> Vec<(String, &Value)> {
      let mut nodes = vec![];
      for (name, schema) in doc["components"]["schemas"].as_object().into_iter().flatten() {
        collect(schema, &format!("{COMPONENTS_PREFIX}{name}"), &mut nodes);
      }
      for (pointer, schema) in operation_schemas(doc) {
        collect(schema, &pointer, &mut nodes);
      }
      nodes
    }

    pub(super) fn open_schemas(doc: &Value) -> Vec<String> {
      schema_nodes(doc)
        .into_iter()
        .filter(|(_, node)| is_open(node))
        .map(|(pointer, _)| pointer)
        .collect()
    }

    pub(super) fn unreachable_components(doc: &Value, roots: &[&str]) -> Vec<String> {
      let schemas = doc["components"]["schemas"].as_object().cloned().unwrap_or_default();
      let mut pending = operation_schemas(doc)
        .into_iter()
        .flat_map(|(_, schema)| references(schema))
        .chain(roots.iter().map(|root| (*root).to_owned()))
        .collect_vec();
      let mut seen = BTreeSet::new();
      while let Some(name) = pending.pop() {
        if seen.insert(name.clone())
          && let Some(schema) = schemas.get(&name)
        {
          pending.extend(references(schema));
        }
      }
      schemas.keys().filter(|name| !seen.contains(*name)).cloned().collect()
    }

    pub(super) fn keyword_locations_in(doc: &Value, keyword: &str, value: &Value) -> Vec<String> {
      schema_nodes(doc)
        .into_iter()
        .filter(|(_, node)| node.get(keyword) == Some(value))
        .map(|(pointer, _)| format!("{pointer}/{keyword}"))
        .collect()
    }

    pub(super) fn enum_spelling_problems(doc: &Value) -> Vec<String> {
      let enums = schema_nodes(doc)
        .into_iter()
        .filter_map(|(pointer, node)| string_enum(node).map(|values| (pointer, values)))
        .collect_vec();
      let mut by_normalized: BTreeMap<Vec<String>, BTreeMap<Vec<String>, String>> = BTreeMap::new();
      for (pointer, values) in &enums {
        let normalized = values.iter().map(|value| normalize(value)).sorted().collect_vec();
        let spelled = values.iter().cloned().sorted().collect_vec();
        by_normalized
          .entry(normalized)
          .or_default()
          .entry(spelled)
          .or_insert_with(|| pointer.clone());
      }
      let differing = by_normalized
        .values()
        .filter(|spellings| spellings.len() > 1)
        .map(|spellings| {
          format!(
            "{} spell the same values differently",
            spellings.values().sorted().join(" and ")
          )
        });
      let mixed = enums
        .iter()
        .filter(|(_, values)| mixes_case_styles(values))
        .map(|(pointer, _)| format!("{pointer} mixes case styles"));
      differing.chain(mixed).sorted().collect()
    }

    fn operation_schemas(doc: &Value) -> Vec<(String, &Value)> {
      let mut schemas = vec![];
      for (path, item) in doc["paths"].as_object().into_iter().flatten() {
        let path = escape(path);
        for method in METHODS {
          let Some(operation) = item.get(*method) else {
            continue;
          };
          let base = format!("#/paths/{path}/{method}");
          for (index, parameter) in operation["parameters"].as_array().into_iter().flatten().enumerate() {
            if let Some(schema) = parameter.get("schema") {
              schemas.push((format!("{base}/parameters/{index}/schema"), schema));
            }
          }
          for (media, content) in operation["requestBody"]["content"].as_object().into_iter().flatten() {
            if let Some(schema) = content.get("schema") {
              schemas.push((format!("{base}/requestBody/content/{}/schema", escape(media)), schema));
            }
          }
          for (status, response) in operation["responses"].as_object().into_iter().flatten() {
            for (media, content) in response["content"].as_object().into_iter().flatten() {
              if let Some(schema) = content.get("schema") {
                schemas.push((
                  format!("{base}/responses/{status}/content/{}/schema", escape(media)),
                  schema,
                ));
              }
            }
          }
        }
      }
      schemas
    }

    fn collect<'a>(schema: &'a Value, pointer: &str, nodes: &mut Vec<(String, &'a Value)>) {
      nodes.push((pointer.to_owned(), schema));
      let Value::Object(object) = schema else {
        return;
      };
      for key in SCHEMA_MAPS {
        for (name, child) in object.get(*key).and_then(Value::as_object).into_iter().flatten() {
          collect(child, &format!("{pointer}/{key}/{}", escape(name)), nodes);
        }
      }
      for key in SCHEMA_LISTS {
        for (index, child) in object
          .get(*key)
          .and_then(Value::as_array)
          .into_iter()
          .flatten()
          .enumerate()
        {
          collect(child, &format!("{pointer}/{key}/{index}"), nodes);
        }
      }
      for key in SCHEMA_SINGLES {
        if let Some(child) = object.get(*key).filter(|child| child.is_object()) {
          collect(child, &format!("{pointer}/{key}"), nodes);
        }
      }
    }

    fn is_open(node: &Value) -> bool {
      match node {
        Value::Bool(open) => *open,
        Value::Object(object) => {
          object.get("additionalProperties") == Some(&Value::Bool(true))
            || object.keys().all(|key| ANNOTATIONS.contains(&key.as_str()))
        },
        _ => false,
      }
    }

    fn references(value: &Value) -> Vec<String> {
      match value {
        Value::Object(object) => object
          .iter()
          .flat_map(|(key, child)| match child {
            Value::String(reference) if key == "$ref" => reference
              .strip_prefix(COMPONENTS_PREFIX)
              .map(ToOwned::to_owned)
              .into_iter()
              .collect_vec(),
            _ => references(child),
          })
          .collect(),
        Value::Array(items) => items.iter().flat_map(references).collect(),
        _ => vec![],
      }
    }

    fn string_enum(node: &Value) -> Option<Vec<String>> {
      if let Some(values) = node.get("enum").and_then(Value::as_array) {
        let strings = values
          .iter()
          .filter_map(Value::as_str)
          .map(ToOwned::to_owned)
          .collect_vec();
        return (strings.len() > 1).then_some(strings);
      }
      let branches = node.get("oneOf").and_then(Value::as_array)?;
      let constants = branches
        .iter()
        .map(|branch| branch.get("const").and_then(Value::as_str).map(ToOwned::to_owned))
        .collect::<Option<Vec<_>>>()?;
      (constants.len() > 1).then_some(constants)
    }

    fn normalize(value: &str) -> String {
      value.to_lowercase().replace(['-', '_'], "")
    }

    fn mixes_case_styles(values: &[String]) -> bool {
      let styles = values.iter().map(|value| case_style(value)).collect::<BTreeSet<_>>();
      let separated = styles.iter().filter(|style| **style != CaseStyle::Lower).collect_vec();
      match separated.as_slice() {
        [] => false,
        [CaseStyle::Kebab | CaseStyle::Snake] => false,
        [_] => styles.contains(&CaseStyle::Lower),
        _ => true,
      }
    }

    fn case_style(value: &str) -> CaseStyle {
      let lower = value.chars().all(|c| c.is_ascii_lowercase() || c.is_ascii_digit());
      let kebab = value
        .chars()
        .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '-');
      let snake = value
        .chars()
        .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_');
      let starts_upper = value.chars().next().is_some_and(|c| c.is_ascii_uppercase());
      if lower {
        CaseStyle::Lower
      } else if kebab {
        CaseStyle::Kebab
      } else if snake {
        CaseStyle::Snake
      } else if starts_upper && value.chars().all(|c| c.is_ascii_alphanumeric()) {
        CaseStyle::Pascal
      } else {
        CaseStyle::Other
      }
    }

    fn escape(segment: &str) -> String {
      segment.replace('~', "~0").replace('/', "~1")
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
    enum CaseStyle {
      Lower,
      Kebab,
      Snake,
      Pascal,
      Other,
    }
  }
}
