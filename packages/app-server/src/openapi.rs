use app_commands::command::{AppCommand, CheckConfigRequest, CheckConfigResponse, CommandOutcome};
use app_commands::config::schema::draft2020_generator;
use app_commands::job::{JobEvent, TerminalEvent};
use eyre::Report;
use schemars::{JsonSchema, Schema};
use serde_json::{Map, Value, json};
use strum::IntoEnumIterator;
use treetime::progress::LogEvent;
use treetime_schema::{ProgressEvent, VersionInfo};
use treetime_utils::{make_error, make_report};

const DEFS_PREFIX: &str = "#/$defs/";
const COMPONENTS_PREFIX: &str = "#/components/schemas/";

pub(crate) fn config_component(command: AppCommand) -> String {
  let name: &str = command.into();
  let mut chars = name.chars();
  let capitalized = chars
    .next()
    .map(|first| first.to_uppercase().chain(chars).collect::<String>())
    .unwrap_or_default();
  format!("{capitalized}Config")
}

pub(crate) fn add_components(doc: &mut Value) -> Result<(), Report> {
  let mut components = Map::new();
  for command in AppCommand::iter() {
    add_root(&mut components, &config_component(command), command.config_schema())?;
  }
  add_type::<CheckConfigRequest>(&mut components)?;
  add_type::<CheckConfigResponse>(&mut components)?;
  add_type::<CommandOutcome>(&mut components)?;
  add_type::<JobEvent>(&mut components)?;
  add_type::<TerminalEvent>(&mut components)?;
  add_type::<ProgressEvent>(&mut components)?;
  add_type::<LogEvent>(&mut components)?;
  add_type::<VersionInfo>(&mut components)?;

  let schemas = doc
    .pointer_mut("/components/schemas")
    .and_then(Value::as_object_mut)
    .ok_or_else(|| make_report!("the OpenAPI document has no component schemas"))?;
  for (name, schema) in components {
    insert_unique(schemas, &name, schema)?;
  }
  Ok(())
}

pub(crate) fn schema_ref(component: &str) -> Value {
  json!({ "$ref": format!("{COMPONENTS_PREFIX}{component}") })
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
