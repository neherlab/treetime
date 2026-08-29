use crate::cli::diagnostics::checks::{
  config_step_names, config_vars, interpolation_diagnostics, pipeline_schema_diagnostics,
  pipeline_structural_diagnostics, schema_diagnostics,
};
use crate::cli::diagnostics::source::{ConfigSource, render_and_bail};
use eyre::Report;
use schemars::Schema;
use serde_json::Value;

/// Diagnose a pipeline config, rendering a batched caret report and returning an error on any problem.
///
/// Shape errors are reported first and on their own: a document whose skeleton is wrong (missing
/// `steps`, a non-mapping step) would make the schema and interpolation passes report cascading
/// noise. Once the shape is sound, structural schema violations and static interpolation-reference
/// errors are gathered together into one batch.
pub fn check_pipeline(source: &ConfigSource, value: &Value) -> Result<(), Report> {
  let structural = pipeline_structural_diagnostics(value);
  if !structural.is_empty() {
    return render_and_bail(source, "invalid pipeline configuration", structural);
  }

  let mut diags = pipeline_schema_diagnostics(value);
  let vars = config_vars(value);
  let step_names = config_step_names(value);
  diags.extend(interpolation_diagnostics(value, &vars, &step_names));
  render_and_bail(source, "invalid pipeline configuration", diags)
}

/// Diagnose a per-command config against its strict schema, rendering a batched caret report.
///
/// Per-command configs have no interpolation and no step structure, so only schema validation
/// applies; templates are not skipped because a command config carries literal values only.
pub fn check_command_config(source: &ConfigSource, value: &Value, schema: &Schema) -> Result<(), Report> {
  let diags = schema_diagnostics(value, schema, false);
  render_and_bail(source, "invalid configuration", diags)
}
