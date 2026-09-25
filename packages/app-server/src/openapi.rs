use app_commands::bridge::error::ErrorResponse;
use app_commands::bridge::operations::DesktopRequest;
use app_commands::check_config::{CheckConfigRequest, CheckConfigResponse};
use app_commands::check_inputs::{CheckInputsRequest, InputFacts};
use app_commands::command::{AppCommand, CommandOutcome};
use app_commands::config::catalog::{SettingCatalog, setting_catalog};
use app_commands::config::cli_flags::annotated_config_schema;
use app_commands::config::schema::draft2020_generator;
use app_commands::job::{IterationEvent, JobEvent, TerminalEvent};
use app_commands::results::clades::{CladeInRuns, CladeRequest};
use app_commands::results::compare::RunComparison;
use app_commands::results::run_results::RunResults;
use app_commands::run_config::{RunConfigRequest, RunConfigResponse};
use app_commands::runs::events::RunEvent;
use app_commands::runs::files::RunFile;
use app_commands::runs::manager::UploadedInput;
use app_commands::runs::record::{
  CancelRunResponse, CreateRunRequest, RunList, RunRecord, RunSummary, StartRunRequest, UpdateRunRequest,
};
use app_commands::runs::setting_differences::SettingDifference;
use app_datasets::DatasetCatalog;
use eyre::Report;
use schemars::{JsonSchema, Schema};
use serde_json::{Map, Value, json};
use strum::IntoEnumIterator;
use treetime::progress::LogEvent;
use treetime_schema::{ProgressEvent, VersionInfo};
use treetime_utils::{make_error, make_report};

const DEFS_PREFIX: &str = "#/$defs/";
const COMPONENTS_PREFIX: &str = "#/components/schemas/";
const SETTING_CATALOG_KEY: &str = "x-setting-catalog";

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
    add_root(
      &mut components,
      &config_component(command),
      annotated_config_schema(command)?,
    )?;
  }
  add_type::<CheckConfigRequest>(&mut components)?;
  add_type::<CheckConfigResponse>(&mut components)?;
  add_type::<RunConfigRequest>(&mut components)?;
  add_type::<RunConfigResponse>(&mut components)?;
  add_type::<CommandOutcome>(&mut components)?;
  add_type::<JobEvent>(&mut components)?;
  add_type::<TerminalEvent>(&mut components)?;
  add_type::<ProgressEvent>(&mut components)?;
  add_type::<LogEvent>(&mut components)?;
  add_type::<VersionInfo>(&mut components)?;
  add_type::<IterationEvent>(&mut components)?;
  add_type::<CheckInputsRequest>(&mut components)?;
  add_type::<InputFacts>(&mut components)?;
  add_type::<DatasetCatalog>(&mut components)?;
  add_type::<CreateRunRequest>(&mut components)?;
  add_type::<StartRunRequest>(&mut components)?;
  add_type::<UpdateRunRequest>(&mut components)?;
  add_type::<RunRecord>(&mut components)?;
  add_type::<RunSummary>(&mut components)?;
  add_type::<RunList>(&mut components)?;
  add_type::<RunEvent>(&mut components)?;
  add_type::<RunFile>(&mut components)?;
  add_type::<UploadedInput>(&mut components)?;
  add_type::<SettingCatalog>(&mut components)?;
  add_type::<SettingDifference>(&mut components)?;
  add_type::<RunResults>(&mut components)?;
  add_type::<RunComparison>(&mut components)?;
  add_type::<CladeRequest>(&mut components)?;
  add_type::<CladeInRuns>(&mut components)?;
  add_type::<ErrorResponse>(&mut components)?;
  add_type::<CancelRunResponse>(&mut components)?;
  add_type::<DesktopRequest>(&mut components)?;

  let schemas = doc
    .as_object_mut()
    .ok_or_else(|| make_report!("the OpenAPI document must be a JSON object"))?
    .entry("components")
    .or_insert_with(|| json!({}))
    .as_object_mut()
    .ok_or_else(|| make_report!("the OpenAPI components must be a JSON object"))?
    .entry("schemas")
    .or_insert_with(|| json!({}))
    .as_object_mut()
    .ok_or_else(|| make_report!("the OpenAPI component schemas must be a JSON object"))?;
  for (name, schema) in components {
    insert_unique(schemas, &name, schema)?;
  }
  Ok(())
}

pub(crate) fn add_setting_catalog(doc: &mut Value) -> Result<(), Report> {
  let catalog = serde_json::to_value(setting_catalog()?)?;
  doc
    .as_object_mut()
    .ok_or_else(|| make_report!("the OpenAPI document must be a JSON object"))?
    .insert(SETTING_CATALOG_KEY.to_owned(), catalog);
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
