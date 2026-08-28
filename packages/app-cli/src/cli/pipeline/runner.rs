use crate::cli::pipeline::resolve::{PipelineDoc, ResolvedPipeline, ResolvedStep, resolve_pipeline};
use crate::cli::pipeline::types::PipelineStepCommand;
use crate::cli::rtt_chart::{write_clock_regression_chart_png, write_clock_regression_chart_svg};
use eyre::Report;
use itertools::Itertools;
use serde_json::{Map, Value};
use std::collections::BTreeSet;
use std::path::Path;
use treetime::commands::ancestral::run::run_ancestral_reconstruction;
use treetime::commands::clock::run::run_clock;
use treetime::commands::mugration::run::run_mugration;
use treetime::commands::optimize::run::run_optimize;
use treetime::commands::prune::run::run_prune;
use treetime::commands::timetree::run::run_timetree_estimation;
use treetime::progress::ProgressSink;
use treetime_utils::io::json::json_or_yaml_read_file;
use treetime_utils::make_error;

/// Load a pipeline config file and resolve it into concrete steps.
///
/// Both JSON and YAML are accepted (format chosen by extension). The document is parsed into a
/// permissive value, then resolved: vars, working directory, per-step output directories, and
/// chained `{{ steps.* }}` references are all made concrete before anything runs.
pub fn load_pipeline(config: &Path) -> Result<ResolvedPipeline, Report> {
  let value: Value = json_or_yaml_read_file(config)?;
  let doc = PipelineDoc::from_value(value)?;
  resolve_pipeline(&doc, &process_env())
}

/// The process environment as a template context object of string values.
fn process_env() -> Value {
  let env = std::env::vars()
    .map(|(key, value)| (key, Value::String(value)))
    .collect::<Map<_, _>>();
  Value::Object(env)
}

/// Select the steps to run, in list order, honoring an optional `--steps` subset.
///
/// An unknown selected name is rejected with a did-you-mean over the pipeline's step names.
pub fn select_steps<'a>(
  pipeline: &'a ResolvedPipeline,
  selected: Option<&BTreeSet<String>>,
) -> Result<Vec<&'a ResolvedStep>, Report> {
  let Some(selected) = selected else {
    return Ok(pipeline.steps.iter().collect());
  };

  let names: BTreeSet<&str> = pipeline.steps.iter().map(|step| step.name.as_str()).collect();
  for name in selected {
    if !names.contains(name.as_str()) {
      let candidates: Vec<&str> = names.iter().copied().collect();
      return make_error!(
        "unknown step `{name}` in --steps; {}",
        crate::cli::pipeline::suggest::suggestion_suffix(name, &candidates)
      );
    }
  }

  Ok(
    pipeline
      .steps
      .iter()
      .filter(|step| selected.contains(&step.name))
      .collect(),
  )
}

/// Run a resolved pipeline, one step at a time, stopping at the first failure.
///
/// Steps run sequentially in the current process; the global thread pool set up by the caller is
/// shared across all of them. When a step fails, no later step runs and the error names the failing
/// step, lists the steps that completed with their output directories, and shows how to resume the
/// remainder with `--steps`.
pub fn run_pipeline(
  pipeline: &ResolvedPipeline,
  selected: Option<&BTreeSet<String>>,
  progress: &dyn ProgressSink,
) -> Result<(), Report> {
  let steps = select_steps(pipeline, selected)?;

  let mut completed: Vec<&ResolvedStep> = Vec::new();
  for (position, step) in steps.iter().enumerate() {
    if let Err(err) = run_step(step, progress) {
      let remaining = steps[position..].iter().map(|step| step.name.as_str()).join(",");
      return Err(err.wrap_err(failure_report(step, &completed, &remaining)));
    }
    completed.push(step);
  }

  Ok(())
}

/// Run one step by dispatching to its command runner.
///
/// Mirrors the per-command dispatch in `main`, including writing the clock regression charts for a
/// clock step. The interactive terminal chart is intentionally skipped: a pipeline is non-interactive.
fn run_step(step: &ResolvedStep, progress: &dyn ProgressSink) -> Result<(), Report> {
  match &step.command {
    PipelineStepCommand::Timetree(args) => run_timetree_estimation(args, progress).map(|_| ()),
    PipelineStepCommand::Optimize(args) => run_optimize(args, progress).map(|_| ()),
    PipelineStepCommand::Prune(args) => run_prune(args, progress).map(|_| ()),
    PipelineStepCommand::Ancestral(args) => run_ancestral_reconstruction(args, progress).map(|_| ()),
    PipelineStepCommand::Mugration(args) => run_mugration(args, progress).map(|_| ()),
    PipelineStepCommand::Clock(args) => {
      let result = run_clock(args, progress)?;
      if let Some(outdir) = &args.output.output_all {
        write_clock_regression_chart_svg(
          &result.regression_results,
          &result.clock_model,
          outdir.join("clock.svg"),
        )?;
        write_clock_regression_chart_png(
          &result.regression_results,
          &result.clock_model,
          outdir.join("clock.png"),
        )?;
      }
      Ok(())
    },
  }
}

/// Compose the mid-pipeline failure message: the failing step, completed steps, and how to resume.
fn failure_report(failed: &ResolvedStep, completed: &[&ResolvedStep], remaining: &str) -> String {
  let done = if completed.is_empty() {
    "none".to_owned()
  } else {
    completed
      .iter()
      .map(|step| format!("`{}`{}", step.name, output_dir_hint(step.outputs.output_all.as_deref())))
      .join(", ")
  };
  format!(
    "pipeline step `{}` failed; completed steps: {done}; resume the remaining steps with --steps={remaining}",
    failed.name
  )
}

/// A ` (outputs in <dir>)` suffix when the step has a known output directory.
fn output_dir_hint(output_all: Option<&Path>) -> String {
  output_all.map_or_else(String::new, |dir| format!(" (outputs in {})", dir.display()))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::cli::pipeline::resolve::PipelineDoc;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use treetime_utils::assert_error;

  fn resolve(config: Value) -> ResolvedPipeline {
    let doc = PipelineDoc::from_value(config).unwrap();
    resolve_pipeline(&doc, &json!({})).unwrap()
  }

  fn config() -> Value {
    json!({
      "output_all": "tmp/run",
      "steps": [
        { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv" } },
        { "name": "anc", "ancestral": { "tree": "{{ steps.tt.outputs.nwk }}" } }
      ]
    })
  }

  // With no subset, every step runs in list order.
  #[test]
  fn test_runner_select_all_steps_in_order() {
    let pipeline = resolve(config());
    let selected = select_steps(&pipeline, None).unwrap();
    assert_eq!(
      vec!["tt", "anc"],
      selected.iter().map(|step| step.name.as_str()).collect::<Vec<_>>()
    );
  }

  // A subset runs only the named steps, still in list order.
  #[test]
  fn test_runner_select_subset_keeps_list_order() {
    let pipeline = resolve(config());
    let selected = select_steps(&pipeline, Some(&btreeset! { "anc".to_owned() })).unwrap();
    assert_eq!(
      vec!["anc"],
      selected.iter().map(|step| step.name.as_str()).collect::<Vec<_>>()
    );
  }

  // An unknown selected step name is rejected with a did-you-mean.
  #[test]
  fn test_runner_select_unknown_step_errors() {
    let pipeline = resolve(config());
    let result = select_steps(&pipeline, Some(&btreeset! { "tta".to_owned() }));
    assert_error!(
      result,
      "unknown step `tta` in --steps; did you mean `tt`? Valid values: `anc`, `tt`"
    );
  }
}
