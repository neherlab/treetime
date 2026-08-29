use crate::cli::pipeline::types::Pipeline;
use clap::ValueEnum;
use eyre::Report;
use log::info;
use schemars::generate::SchemaSettings;
use schemars::transform::{Transform, transform_subschemas};
use schemars::{JsonSchema, Schema, SchemaGenerator};
use serde_json::{Value, json};
use std::io::Write;
use std::path::{Path, PathBuf};
use treetime::commands::ancestral::args::TreetimeAncestralArgs;
use treetime::commands::clock::args::TreetimeClockArgs;
use treetime::commands::mugration::args::TreetimeMugrationArgs;
use treetime::commands::optimize::args::TreetimeOptimizeArgs;
use treetime::commands::prune::args::TreetimePruneArgs;
use treetime::commands::timetree::args::TreetimeTimetreeArgs;
use treetime_schema::{TreetimeSchemaFormat, generate_schema as generate_data_schema};
use treetime_utils::io::json::{JsonPretty, json_write_str};

/// A schema the `schema` subcommand can emit.
///
/// Three groups share one selector: the runtime data-type schemas (delegated to `treetime-schema`),
/// the pipeline config schema, and the per-command config schemas. The pipeline and per-command
/// schemas live here rather than in `treetime-schema` because they reference the command argument
/// types, which `treetime-schema` cannot depend on without a crate cycle.
#[derive(Debug, Clone, Copy, Default, ValueEnum, serde::Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum SchemaTarget {
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
  /// The data-type schema this target delegates to `treetime-schema`, if any.
  const fn data_format(self) -> Option<TreetimeSchemaFormat> {
    match self {
      SchemaTarget::VersionInfo => Some(TreetimeSchemaFormat::VersionInfo),
      SchemaTarget::ProgressEvent => Some(TreetimeSchemaFormat::ProgressEvent),
      SchemaTarget::ErrorResponse => Some(TreetimeSchemaFormat::ErrorResponse),
      _ => None,
    }
  }

  /// Default output filename for this target, or `None` for the `all` aggregate.
  const fn default_filename(self) -> Option<&'static str> {
    match self {
      SchemaTarget::All => None,
      SchemaTarget::VersionInfo => Some("version-info.schema.json"),
      SchemaTarget::ProgressEvent => Some("progress-event.schema.json"),
      SchemaTarget::ErrorResponse => Some("error-response.schema.json"),
      SchemaTarget::Pipeline => Some("pipeline.schema.json"),
      SchemaTarget::Timetree => Some("timetree.schema.json"),
      SchemaTarget::Optimize => Some("optimize.schema.json"),
      SchemaTarget::Prune => Some("prune.schema.json"),
      SchemaTarget::Ancestral => Some("ancestral.schema.json"),
      SchemaTarget::Clock => Some("clock.schema.json"),
      SchemaTarget::Mugration => Some("mugration.schema.json"),
    }
  }
}

/// Write the schema (or schemas) selected by `target` to `output`.
///
/// `all` writes every schema, using each one's default filename, into the directory `output` (the
/// current directory when omitted). A single target writes to `output` when given, or to stdout.
pub fn generate_schema(target: SchemaTarget, output: Option<&PathBuf>) -> Result<(), Report> {
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

/// Every concrete (non-aggregate) target, in declared order.
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

/// Emit one concrete target to `output`.
fn generate_one(target: SchemaTarget, output: &Path) -> Result<(), Report> {
  if let Some(format) = target.data_format() {
    let path = output.to_path_buf();
    return generate_data_schema(&format, Some(&path));
  }

  let schema = match target {
    SchemaTarget::Pipeline => pipeline_schema(),
    SchemaTarget::Timetree => command_schema::<TreetimeTimetreeArgs>(),
    SchemaTarget::Optimize => command_schema::<TreetimeOptimizeArgs>(),
    SchemaTarget::Prune => command_schema::<TreetimePruneArgs>(),
    SchemaTarget::Ancestral => command_schema::<TreetimeAncestralArgs>(),
    SchemaTarget::Clock => command_schema::<TreetimeClockArgs>(),
    SchemaTarget::Mugration => command_schema::<TreetimeMugrationArgs>(),
    SchemaTarget::All | SchemaTarget::VersionInfo | SchemaTarget::ProgressEvent | SchemaTarget::ErrorResponse => {
      unreachable!("aggregate and data-type targets are handled earlier")
    },
  };

  write_schema(&schema, output)
}

/// The pipeline config schema, with every scalar leaf loosened to also accept a `{{ ... }}` template.
///
/// A whole-value template renders to a string even where the field is typed (`clock_rate: "{{ vars.rate }}"`),
/// so the editor must accept a template string in place of any scalar. Only the pipeline schema is
/// loosened this way; the per-command schemas stay strict, because interpolation is a pipeline feature.
pub fn pipeline_schema() -> Schema {
  let mut schema = draft2020_generator().into_root_schema_for::<Pipeline>();
  AllowTemplateStrings.transform(&mut schema);
  schema
}

/// A strict per-command config schema (no template loosening).
pub fn command_schema<T: JsonSchema>() -> Schema {
  draft2020_generator().into_root_schema_for::<T>()
}

/// A draft 2020-12 generator, matching the dialect `jsonschema` validates against by default.
fn draft2020_generator() -> SchemaGenerator {
  SchemaSettings::draft2020_12().into_generator()
}

/// Serialize a schema to `output` (a file, or stdout for `-`), creating parent directories as needed.
fn write_schema(schema: &Schema, output: &Path) -> Result<(), Report> {
  let json = json_write_str(schema, JsonPretty(true))?;
  if output == Path::new("-") {
    std::io::stdout().write_all(json.as_bytes())?;
  } else {
    if let Some(parent) = output.parent() {
      std::fs::create_dir_all(parent)?;
    }
    std::fs::write(output, json)?;
    info!("Wrote JSON schema to '{}'", output.display());
  }
  Ok(())
}

/// Rewrites every scalar leaf of a schema into `anyOf: [<original>, <template string>]`.
///
/// The template branch is a string constrained to contain a `{{ ... }}` expression, so an editor
/// accepts a whole-value template in a numeric, boolean, or string field while still validating a
/// literal value against the original type. Subschemas are transformed first so a wrapped leaf is
/// never re-wrapped.
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

  /// The regex a template-string branch carries, as it appears in the schema value (one backslash).
  const TEMPLATE_PATTERN: &str = r"\{\{.*\}\}";

  // The pipeline schema loosens scalar leaves so a whole-value template is accepted where a typed
  // value is expected: the step `name` leaf becomes `anyOf: [string, template string]`.
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

  // A per-command schema stays strict: interpolation is a pipeline feature, so no scalar leaf is
  // loosened and no template pattern appears anywhere in the document.
  #[test]
  fn test_schema_command_is_strict_without_templates() {
    let schema = serde_json::to_value(command_schema::<TreetimeAncestralArgs>()).unwrap();
    assert!(
      !helpers::contains_template_pattern(&schema),
      "a per-command schema must not loosen any leaf to a template string"
    );
  }

  // The pipeline schema does apply the loosening (guards the negative test above against a transform
  // that silently stopped firing).
  #[test]
  fn test_schema_pipeline_contains_template_pattern() {
    let schema = serde_json::to_value(pipeline_schema()).unwrap();
    assert!(helpers::contains_template_pattern(&schema));
  }

  mod helpers {
    use serde_json::Value;

    /// Whether any node in the schema carries the template-string `pattern`.
    pub fn contains_template_pattern(value: &Value) -> bool {
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
