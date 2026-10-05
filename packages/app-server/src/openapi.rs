use aide::openapi::{Components, OpenApi, SchemaObject};
use app_commands::command::AppCommand;
use app_commands::config::catalog::{SettingCatalog, setting_catalog};
use app_commands::config::cli_flags::annotated_config_schema;
use app_commands::config::schema::{SCHEMA_KEY, draft2020_generator};
use eyre::Report;
use heck::ToUpperCamelCase;
use indexmap::IndexMap;
use itertools::{Itertools, izip};
use schemars::{JsonSchema, Schema};
use serde_json::{Map, Value, json};
use std::collections::BTreeSet;
use strum::IntoEnumIterator;
use treetime_utils::make_error;

const DEFS_PREFIX: &str = "#/$defs/";
const COMPONENTS_PREFIX: &str = "#/components/schemas/";
const SETTING_CATALOG_KEY: &str = "x-setting-catalog";

const UNION_ROOT_KEYS: &[&str] = &["description", "oneOf", "type", "properties", "required"];

const SCHEMA_MAPS: &[&str] = &["properties", "$defs", "patternProperties"];

const SCHEMA_LISTS: &[&str] = &["anyOf", "oneOf", "allOf", "prefixItems"];

const SCHEMA_SINGLES: &[&str] = &["items", "additionalProperties", "not", "if", "then", "else"];

pub(crate) fn add_components(api: &mut OpenApi) -> Result<(), Report> {
  let mut components = Map::new();
  add_type::<SettingCatalog>(&mut components)?;

  let schemas = &mut api.components.get_or_insert_with(Components::default).schemas;
  replace_config_components(schemas, &mut components)?;
  for (name, schema) in components {
    let schema = Schema::try_from(schema)?;
    match schemas.get(&name) {
      Some(existing) if existing.json_schema != schema => {
        return make_error!("two different schemas are named `{name}`; rename one of the types");
      },
      Some(_) => {},
      None => {
        schemas.insert(
          name,
          SchemaObject {
            json_schema: schema,
            example: None,
            external_docs: None,
          },
        );
      },
    }
  }
  Ok(())
}

pub(crate) fn add_setting_catalog(api: &mut OpenApi) -> Result<(), Report> {
  api.extensions.insert(
    SETTING_CATALOG_KEY.to_owned(),
    serde_json::to_value(setting_catalog()?)?,
  );
  Ok(())
}

pub(crate) fn add_discriminators(api: &mut OpenApi) -> Result<(), Report> {
  let schemas = &mut api.components.get_or_insert_with(Components::default).schemas;
  let unions = schemas
    .iter()
    .filter_map(|(name, schema)| tagged_union(name, &schema.json_schema).transpose())
    .collect::<Result<Vec<_>, _>>()?;
  for union in unions {
    for (variant_name, variant) in union.variants {
      if schemas.contains_key(&variant_name) {
        return make_error!(
          "the variant `{variant_name}` of the tagged union `{}` has the name of another schema; rename one of the types",
          union.name
        );
      }
      schemas.insert(
        variant_name,
        SchemaObject {
          json_schema: variant,
          example: None,
          external_docs: None,
        },
      );
    }
    schemas[&union.name].json_schema = union.schema;
  }
  Ok(())
}

struct TaggedUnion {
  name: String,
  schema: Schema,
  variants: Vec<(String, Schema)>,
}

fn tagged_union(name: &str, schema: &Schema) -> Result<Option<TaggedUnion>, Report> {
  let Some(root) = schema.as_object() else {
    return Ok(None);
  };
  let Some(Value::Array(one_of)) = root.get("oneOf") else {
    return Ok(None);
  };
  let Some(variants) = one_of
    .iter()
    .map(|variant| {
      variant
        .as_object()
        .filter(|variant| variant.get("type") == Some(&json!("object")))
    })
    .collect::<Option<Vec<_>>>()
  else {
    return Ok(None);
  };
  let Some((tag, values)) = tag_property(name, &variants)? else {
    return Ok(None);
  };
  if let Some((key, _)) = root.iter().find(|(key, value)| {
    !UNION_ROOT_KEYS.contains(&key.as_str()) || (key.as_str() == "type" && *value != &json!("object"))
  }) {
    return make_error!(
      "the tagged union `{name}` has the keyword `{key}` next to `oneOf`, which its variants cannot take over"
    );
  }
  let empty = Map::new();
  let shared_properties = root.get("properties").and_then(Value::as_object).unwrap_or(&empty);
  let shared_required = root
    .get("required")
    .and_then(Value::as_array)
    .map_or(&[][..], Vec::as_slice);

  let mut refs = vec![];
  let mut mapping = Map::new();
  let mut named_variants = vec![];
  for (variant, value) in izip!(variants, values) {
    let variant_name = format!("{name}{}", value.to_upper_camel_case());
    let reference = format!("{COMPONENTS_PREFIX}{variant_name}");
    refs.push(json!({ "$ref": reference }));
    mapping.insert(value.to_owned(), json!(reference));
    let variant = with_shared_fields(name, variant, shared_properties, shared_required)?;
    named_variants.push((variant_name, Schema::try_from(Value::Object(variant))?));
  }

  let mut union = Map::new();
  if let Some(description) = root.get("description") {
    union.insert("description".to_owned(), description.clone());
  }
  union.insert("oneOf".to_owned(), Value::Array(refs));
  union.insert(
    "discriminator".to_owned(),
    json!({ "propertyName": tag, "mapping": mapping }),
  );
  Ok(Some(TaggedUnion {
    name: name.to_owned(),
    schema: Schema::try_from(Value::Object(union))?,
    variants: named_variants,
  }))
}

fn tag_property<'a>(name: &str, variants: &[&'a Map<String, Value>]) -> Result<Option<(String, Vec<&'a str>)>, Report> {
  let Some(first) = variants
    .first()
    .and_then(|variant| variant.get("properties"))
    .and_then(Value::as_object)
  else {
    return Ok(None);
  };
  let tags = first
    .keys()
    .filter(|property| {
      variants.iter().all(|variant| {
        tag_value(variant, property).is_some()
          && variant
            .get("required")
            .and_then(Value::as_array)
            .is_some_and(|required| required.contains(&json!(property)))
      })
    })
    .collect_vec();
  match tags.as_slice() {
    [] => Ok(None),
    [tag] => {
      let values = variants
        .iter()
        .filter_map(|variant| tag_value(variant, tag))
        .collect_vec();
      if values.iter().all_unique() {
        Ok(Some(((*tag).clone(), values)))
      } else {
        make_error!("the tagged union `{name}` has two variants with the same `{tag}`")
      }
    },
    _ => make_error!(
      "the union `{name}` can be told apart by each of {}; the discriminator needs exactly one tag property",
      tags.iter().map(|tag| format!("`{tag}`")).join(", ")
    ),
  }
}

fn tag_value<'a>(variant: &'a Map<String, Value>, tag: &str) -> Option<&'a str> {
  variant.get("properties")?.get(tag)?.get("const")?.as_str()
}

fn with_shared_fields(
  name: &str,
  variant: &Map<String, Value>,
  shared_properties: &Map<String, Value>,
  shared_required: &[Value],
) -> Result<Map<String, Value>, Report> {
  let mut variant = variant.clone();
  if shared_properties.is_empty() && shared_required.is_empty() {
    return Ok(variant);
  }
  if variant.get("additionalProperties") == Some(&Value::Bool(false)) {
    return make_error!("the tagged union `{name}` has fields next to a variant that allows no other fields");
  }
  let mut properties = shared_properties.clone();
  if let Some(own) = variant.get("properties").and_then(Value::as_object) {
    for (property, schema) in own {
      match properties.get(property) {
        Some(shared) if shared != schema => {
          return make_error!("the tagged union `{name}` describes the field `{property}` twice, differently");
        },
        _ => {
          properties.insert(property.clone(), schema.clone());
        },
      }
    }
  }
  let own_required = variant
    .get("required")
    .and_then(Value::as_array)
    .cloned()
    .unwrap_or_default();
  let required = shared_required
    .iter()
    .chain(&own_required)
    .filter_map(Value::as_str)
    .unique()
    .map(|property| json!(property))
    .collect_vec();
  variant.insert("properties".to_owned(), Value::Object(properties));
  variant.insert("required".to_owned(), Value::Array(required));
  Ok(variant)
}

fn replace_config_components(
  schemas: &mut IndexMap<String, SchemaObject>,
  components: &mut Map<String, Value>,
) -> Result<(), Report> {
  let configs = AppCommand::iter()
    .map(|command| (command, command.config_schema_name()))
    .filter(|(_, name)| schemas.contains_key(name))
    .collect_vec();
  let config_names = configs.iter().map(|(_, name)| name.clone()).collect::<BTreeSet<_>>();
  let below_configs = reachable(schemas, &config_names, &BTreeSet::new());
  let others = schemas
    .keys()
    .filter(|name| !below_configs.contains(*name))
    .cloned()
    .collect::<BTreeSet<_>>();
  let below_others = reachable(schemas, &others, &config_names);
  for name in below_configs.difference(&below_others) {
    schemas.shift_remove(name);
  }
  for (command, name) in configs {
    let mut schema = annotated_config_schema(command)?.to_value();
    if let Some(Value::Object(properties)) = schema.get_mut("properties") {
      properties.remove(SCHEMA_KEY);
    }
    strip_defaults(&mut schema);
    add_root(components, &name, Schema::try_from(schema)?)?;
  }
  Ok(())
}

fn reachable(
  schemas: &IndexMap<String, SchemaObject>,
  roots: &BTreeSet<String>,
  stop: &BTreeSet<String>,
) -> BTreeSet<String> {
  let mut seen = BTreeSet::new();
  let mut pending = roots.iter().cloned().collect_vec();
  while let Some(name) = pending.pop() {
    if !seen.insert(name.clone()) {
      continue;
    }
    if let Some(schema) = schemas.get(&name) {
      let mut refs = vec![];
      component_refs(schema.json_schema.as_value(), &mut refs);
      pending.extend(refs.into_iter().filter(|reference| !stop.contains(reference)));
    }
  }
  seen
}

fn component_refs(value: &Value, refs: &mut Vec<String>) {
  match value {
    Value::Object(object) => {
      for (key, child) in object {
        match child {
          Value::String(reference) if key == "$ref" => {
            if let Some(name) = reference.strip_prefix(COMPONENTS_PREFIX) {
              refs.push(name.to_owned());
            }
          },
          _ => component_refs(child, refs),
        }
      }
    },
    Value::Array(items) => items.iter().for_each(|item| component_refs(item, refs)),
    _ => {},
  }
}

fn strip_defaults(schema: &mut Value) {
  let Value::Object(object) = schema else {
    return;
  };
  object.remove("default");
  for key in SCHEMA_MAPS {
    if let Some(Value::Object(children)) = object.get_mut(*key) {
      children.values_mut().for_each(strip_defaults);
    }
  }
  for key in SCHEMA_LISTS {
    if let Some(Value::Array(children)) = object.get_mut(*key) {
      children.iter_mut().for_each(strip_defaults);
    }
  }
  for key in SCHEMA_SINGLES {
    if let Some(child) = object.get_mut(*key) {
      strip_defaults(child);
    }
  }
}

fn add_type<T: JsonSchema>(components: &mut Map<String, Value>) -> Result<(), Report> {
  let schema = draft2020_generator().into_root_schema_for::<T>();
  add_root(components, &T::schema_name(), schema)
}

fn add_root(components: &mut Map<String, Value>, name: &str, schema: Schema) -> Result<(), Report> {
  let mut root = schema.to_value();
  let defs = match &mut root {
    Value::Object(object) => {
      object.remove("$schema");
      object.remove("title");
      object.remove("$defs")
    },
    _ => return make_error!("the schema of `{name}` is not an object"),
  };
  rewrite_refs(&mut root);
  if let Some(Value::Object(defs)) = defs {
    for (def_name, mut def) in defs {
      rewrite_refs(&mut def);
      insert_unique(components, &def_name, def)?;
    }
  }
  insert_unique(components, name, root)
}

fn insert_unique(components: &mut Map<String, Value>, name: &str, schema: Value) -> Result<(), Report> {
  match components.get(name) {
    Some(existing) if *existing != schema => {
      make_error!("two different schemas are named `{name}`; rename one of the types")
    },
    Some(_) => Ok(()),
    None => {
      components.insert(name.to_owned(), schema);
      Ok(())
    },
  }
}

fn rewrite_refs(value: &mut Value) {
  match value {
    Value::Object(object) => {
      for (key, child) in object.iter_mut() {
        match child {
          Value::String(reference) if key == "$ref" => {
            if let Some(def) = reference.strip_prefix(DEFS_PREFIX) {
              *reference = format!("{COMPONENTS_PREFIX}{def}");
            }
          },
          _ => rewrite_refs(child),
        }
      }
    },
    Value::Array(items) => items.iter_mut().for_each(rewrite_refs),
    _ => {},
  }
}
