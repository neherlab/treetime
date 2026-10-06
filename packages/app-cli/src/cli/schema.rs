use crate::cli::pipeline::types::Pipeline;
use app_commands::command::AppCommand;
use app_commands::config::cli_flags::annotated_config_schema;
use app_commands::config::schema::draft2020_generator;
use clap::ValueEnum;
use deser::Serialize;
use deser::adapters::As;
use deser_serde::Serde;
use eyre::Report;
use log::info;
use schemars::Schema;
use schemars::transform::{Transform, transform_subschemas};
use serde_json::{Value, json};
use std::path::{Path, PathBuf};
use treetime_schema::TreetimeSchemaFormat;
use treetime_utils::io::file::is_path_stdout;
use treetime_utils::io::json::{JsonPretty, json_write_file};

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
pub(crate) fn generate_schema(target: SchemaTarget, output: Option<&PathBuf>) -> Result<(), Report> {
  if matches!(target, SchemaTarget::All) {
    let dir = output.map_or_else(|| PathBuf::from("."), Clone::clone);
    for one in all_targets() {
      let filename = one.default_filename().expect("non-aggregate target has a filename");
      generate_one(one, &dir.join(filename))?;
    }
    return Ok(());
  }

  let path = output.map_or_else(|| PathBuf::from("-"), Clone::clone);
  generate_one(target, &path)
}

fn all_targets() -> impl Iterator<Item = SchemaTarget> {
  [
    SchemaTarget::VersionInfo,
    SchemaTarget::ProgressEvent,
    SchemaTarget::ErrorResponse,
    SchemaTarget::Pipeline,
    SchemaTarget::Timetree,
    SchemaTarget::Optimize,
    SchemaTarget::Prune,
    SchemaTarget::Ancestral,
    SchemaTarget::Homoplasy,
    SchemaTarget::Clock,
    SchemaTarget::Mugration,
  ]
  .into_iter()
}

fn generate_one(target: SchemaTarget, output: &Path) -> Result<(), Report> {
  if let Some(format) = target.data_format() {
    let path = output.to_path_buf();
    return treetime_schema::generate_schema(&format, Some(&path));
  }

  let schema = match target {
    SchemaTarget::Pipeline => pipeline_schema(),
    SchemaTarget::Timetree => annotated_config_schema(AppCommand::Timetree)?,
    SchemaTarget::Optimize => annotated_config_schema(AppCommand::Optimize)?,
    SchemaTarget::Prune => annotated_config_schema(AppCommand::Prune)?,
    SchemaTarget::Ancestral => annotated_config_schema(AppCommand::Ancestral)?,
    SchemaTarget::Homoplasy => annotated_config_schema(AppCommand::Homoplasy)?,
    SchemaTarget::Clock => annotated_config_schema(AppCommand::Clock)?,
    SchemaTarget::Mugration => annotated_config_schema(AppCommand::Mugration)?,
    SchemaTarget::All | SchemaTarget::VersionInfo | SchemaTarget::ProgressEvent | SchemaTarget::ErrorResponse => {
      unreachable!("aggregate and data-type targets are handled earlier")
    },
  };

  write_schema(&schema, output)
}

#[derive(Debug, Clone, Copy, Default, ValueEnum, Serialize)]
#[deser(rename_all = "kebab-case")]
pub(crate) enum SchemaTarget {
  #[default]
  All,
  VersionInfo,
  ProgressEvent,
  ErrorResponse,
  Pipeline,
  Timetree,
  Optimize,
  Prune,
  Ancestral,
  Homoplasy,
  Clock,
  Mugration,
}

impl SchemaTarget {
  const fn data_format(self) -> Option<TreetimeSchemaFormat> {
    match self {
      SchemaTarget::VersionInfo => Some(TreetimeSchemaFormat::VersionInfo),
      SchemaTarget::ProgressEvent => Some(TreetimeSchemaFormat::ProgressEvent),
      SchemaTarget::ErrorResponse => Some(TreetimeSchemaFormat::ErrorResponse),
      _ => None,
    }
  }

  const fn default_filename(self) -> Option<&'static str> {
    match self {
      SchemaTarget::All => None,
      SchemaTarget::VersionInfo => Some("version-info.schema.json"),
      SchemaTarget::ProgressEvent => Some("progress-event.schema.json"),
      SchemaTarget::ErrorResponse => Some("error-response.schema.json"),
      SchemaTarget::Pipeline => Some("input-config-pipeline.schema.json"),
      SchemaTarget::Timetree => Some("input-config-timetree.schema.json"),
      SchemaTarget::Optimize => Some("input-config-optimize.schema.json"),
      SchemaTarget::Prune => Some("input-config-prune.schema.json"),
      SchemaTarget::Ancestral => Some("input-config-ancestral.schema.json"),
      SchemaTarget::Homoplasy => Some("input-config-homoplasy.schema.json"),
      SchemaTarget::Clock => Some("input-config-clock.schema.json"),
      SchemaTarget::Mugration => Some("input-config-mugration.schema.json"),
    }
  }
}

fn pipeline_schema() -> Schema {
  let mut schema = draft2020_generator().into_root_schema_for::<Pipeline>();
  AllowTemplateStrings.transform(&mut schema);
  schema
}

pub(crate) fn command_schema_for(tag: &str) -> Option<Schema> {
  tag.parse::<AppCommand>().ok().map(AppCommand::config_schema)
}

fn write_schema(schema: &Schema, output: &Path) -> Result<(), Report> {
  json_write_file(output, &As::<_, Serde>::new(schema), JsonPretty(true))?;
  if !is_path_stdout(output) {
    info!("Wrote JSON schema to '{}'", output.display());
  }
  Ok(())
}

struct AllowTemplateStrings;

impl Transform for AllowTemplateStrings {
  fn transform(&mut self, schema: &mut Schema) {
    transform_subschemas(self, schema);

    let is_scalar_leaf = schema
      .get("type")
      .and_then(Value::as_str)
      .is_some_and(|ty| matches!(ty, "string" | "number" | "integer" | "boolean"));
    if !is_scalar_leaf {
      return;
    }

    let object = schema.ensure_object();
    let original = Value::Object(object.clone());
    object.clear();
    object.insert(
      "anyOf".to_owned(),
      json!([original, { "type": "string", "pattern": "\\{\\{.*\\}\\}" }]),
    );
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use std::fs;
  use tempfile::tempdir;

  const TEMPLATE_PATTERN: &str = r"\{\{.*\}\}";

  #[test]
  fn test_schema_committed_files_match_generated() {
    let dir = tempdir().unwrap();
    generate_schema(SchemaTarget::All, Some(&dir.path().to_path_buf())).unwrap();

    let committed_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../schemas");
    for target in all_targets() {
      let filename = target.default_filename().expect("non-aggregate target has a filename");
      let generated: Value = serde_json::from_str(&fs::read_to_string(dir.path().join(filename)).unwrap()).unwrap();
      let committed: Value = serde_json::from_str(&fs::read_to_string(committed_dir.join(filename)).unwrap()).unwrap();
      assert_eq!(
        committed, generated,
        "committed schema `{filename}` is stale; regenerate with `treetime schema --for all -o packages/schemas`"
      );
    }
  }

  #[test]
  fn test_schema_no_property_defaults_to_null() {
    let dir = tempdir().unwrap();
    generate_schema(SchemaTarget::All, Some(&dir.path().to_path_buf())).unwrap();
    let null_defaults = all_targets()
      .flat_map(|target| {
        let filename = target.default_filename().expect("non-aggregate target has a filename");
        let schema: Value = serde_json::from_str(&fs::read_to_string(dir.path().join(filename)).unwrap()).unwrap();
        helpers::null_defaults(&schema, filename)
      })
      .collect::<Vec<_>>();
    assert_eq!(Vec::<String>::new(), null_defaults);
  }

  #[test]
  fn test_schema_pipeline_loosens_scalar_leaf_to_template() {
    let schema = serde_json::to_value(pipeline_schema()).unwrap();
    let branches = &schema["$defs"]["PipelineStep"]["properties"]["name"]["anyOf"];
    let patterns: Vec<&str> = branches
      .as_array()
      .expect("name leaf is an anyOf")
      .iter()
      .filter_map(|branch| branch.get("pattern").and_then(Value::as_str))
      .collect();
    assert_eq!(vec![TEMPLATE_PATTERN], patterns);
  }

  #[test]
  fn test_schema_command_is_strict_without_templates() {
    let schema = serde_json::to_value(AppCommand::Ancestral.config_schema()).unwrap();
    assert!(
      !helpers::contains_template_pattern(&schema),
      "a per-command schema must not loosen any leaf to a template string"
    );
  }

  #[test]
  fn test_schema_pipeline_contains_template_pattern() {
    let schema = serde_json::to_value(pipeline_schema()).unwrap();
    assert!(helpers::contains_template_pattern(&schema));
  }

  mod helpers {
    use serde_json::Value;

    pub(super) fn contains_template_pattern(value: &Value) -> bool {
      match value {
        Value::Object(map) => {
          let here = map
            .get("pattern")
            .and_then(Value::as_str)
            .is_some_and(|pattern| pattern == super::TEMPLATE_PATTERN);
          here || map.values().any(contains_template_pattern)
        },
        Value::Array(items) => items.iter().any(contains_template_pattern),
        _ => false,
      }
    }

    pub(super) fn null_defaults(value: &Value, pointer: &str) -> Vec<String> {
      match value {
        Value::Object(object) => object
          .iter()
          .flat_map(|(key, child)| {
            let here = format!("{pointer}/{key}");
            if key == "default" && child.is_null() {
              vec![here]
            } else {
              null_defaults(child, &here)
            }
          })
          .collect(),
        Value::Array(items) => items
          .iter()
          .enumerate()
          .flat_map(|(index, item)| null_defaults(item, &format!("{pointer}/{index}")))
          .collect(),
        _ => vec![],
      }
    }
  }
}
