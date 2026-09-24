use crate::cli::diagnostics::checks::{
  config_step_names, config_vars, interpolation_diagnostics, pipeline_schema_diagnostics,
  pipeline_structural_diagnostics,
};
use app_commands::config::source::{ConfigSource, render_and_bail};
use eyre::Report;
use serde_json::Value;

pub(crate) fn check_pipeline(source: &ConfigSource, value: &Value) -> Result<(), Report> {
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
