use crate::cli::pipeline::types::{Pipeline, SCHEMA_KEY};
use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
use crate::commands::clock::args::TreetimeClockArgsRaw;
use crate::commands::mugration::args::TreetimeMugrationArgsRaw;
use crate::commands::optimize::args::TreetimeOptimizeArgsRaw;
use crate::commands::prune::args::TreetimePruneArgsRaw;
use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
use clap::ValueEnum;
use eyre::{Report, WrapErr};
use log::info;
use schemars::generate::SchemaSettings;
use schemars::transform::{Transform, transform_subschemas};
use schemars::{JsonSchema, Schema, SchemaGenerator};
use serde_json::{Value, json};
use std::io::Write;
use std::path::{Path, PathBuf};
use treetime_schema::TreetimeSchemaFormat;
use treetime_utils::io::json::{JsonPretty, json_write_str};

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
    SchemaTarget::Timetree => command_schema::<TreetimeTimetreeArgsRaw>(),
    SchemaTarget::Optimize => command_schema::<TreetimeOptimizeArgsRaw>(),
    SchemaTarget::Prune => command_schema::<TreetimePruneArgsRaw>(),
    SchemaTarget::Ancestral => command_schema::<TreetimeAncestralArgsRaw>(),
    SchemaTarget::Clock => command_schema::<TreetimeClockArgsRaw>(),
    SchemaTarget::Mugration => command_schema::<TreetimeMugrationArgsRaw>(),
    SchemaTarget::All | SchemaTarget::VersionInfo | SchemaTarget::ProgressEvent | SchemaTarget::ErrorResponse => {
      unreachable!("aggregate and data-type targets are handled earlier")
    },
  };

  write_schema(&schema, output)
}

#[derive(Debug, Clone, Copy, Default, ValueEnum, serde::Serialize)]
#[serde(rename_all = "kebab-case")]
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
  Some(match tag {
    "timetree" => command_schema::<TreetimeTimetreeArgsRaw>(),
    "optimize" => command_schema::<TreetimeOptimizeArgsRaw>(),
    "prune" => command_schema::<TreetimePruneArgsRaw>(),
    "ancestral" => command_schema::<TreetimeAncestralArgsRaw>(),
    "clock" => command_schema::<TreetimeClockArgsRaw>(),
    "mugration" => command_schema::<TreetimeMugrationArgsRaw>(),
    _ => return None,
  })
}

pub(crate) fn command_schema<T: JsonSchema>() -> Schema {
  let mut schema = draft2020_generator().into_root_schema_for::<T>();
  allow_schema_ref(&mut schema);
  schema
}

fn allow_schema_ref(schema: &mut Schema) {
  let object = schema.ensure_object();
  let properties = object.entry("properties").or_insert_with(|| json!({}));
  if let Some(properties) = properties.as_object_mut() {
    properties.insert(
      SCHEMA_KEY.to_owned(),
      json!({
        "type": "string",
        "description": "Path or URL of the JSON schema for this config; used by editors and ignored by the loader."
      }),
    );
  }
}

fn draft2020_generator() -> SchemaGenerator {
  SchemaSettings::draft2020_12().into_generator()
}

fn write_schema(schema: &Schema, output: &Path) -> Result<(), Report> {
  let json = json_write_str(schema, JsonPretty(true))?;
  if output == Path::new("-") {
    std::io::stdout()
      .write_all(json.as_bytes())
      .wrap_err("When writing the JSON schema to standard output")?;
  } else {
    if let Some(parent) = output.parent() {
      std::fs::create_dir_all(parent).wrap_err_with(|| format!("When creating directory '{}'", parent.display()))?;
    }
    std::fs::write(output, json).wrap_err_with(|| format!("When writing JSON schema file '{}'", output.display()))?;
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
    let schema = serde_json::to_value(command_schema::<TreetimeAncestralArgsRaw>()).unwrap();
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
  }
}
