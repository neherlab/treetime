use crate::cli::diagnostics::entry::check_pipeline;
use crate::cli::diagnostics::source::{ConfigSource, parse_config_document};
use crate::cli::pipeline::check::print_pipeline_plan;
use crate::cli::pipeline::resolve::{PipelineDoc, ResolvedPipeline, ResolvedStep, resolve_pipeline};
use crate::cli::pipeline::safety::validate_plan;
use crate::cli::pipeline::suggest::suggestion_suffix;
use crate::cli::pipeline::types::PipelineStepCommand;
use crate::cli::treetime_cli::TreetimePipelineArgs;
use crate::commands::ancestral::args::TreetimeAncestralArgs;
use crate::commands::ancestral::run::run_ancestral_reconstruction;
use crate::commands::clock::args::TreetimeClockArgs;
use crate::commands::clock::run::run_clock;
use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::mugration::run::run_mugration;
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::optimize::run::run_optimize;
use crate::commands::prune::args::TreetimePruneArgs;
use crate::commands::prune::run::run_prune;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::run::run_timetree_estimation;
use eyre::Report;
use itertools::Itertools;
use serde_json::{Map, Value};
use std::collections::BTreeSet;
use std::path::Path;
use treetime::cancel::NoopCancel;
use treetime::progress::ProgressSink;
use treetime_utils::io::fs::read_file_to_string;
use treetime_utils::make_error;

pub(crate) fn run_pipeline_command(args: &TreetimePipelineArgs, progress: &dyn ProgressSink) -> Result<(), Report> {
  let pipeline = load_pipeline(&args.config)?;
  let selected = (!args.steps.is_empty()).then(|| args.steps.iter().cloned().collect::<BTreeSet<String>>());
  validate_plan(&pipeline, selected.as_ref())?;
  if args.check {
    print_pipeline_plan(&pipeline, selected.as_ref())
  } else {
    run_pipeline(&pipeline, selected.as_ref(), progress)
  }
}

pub(crate) fn load_pipeline(config: &Path) -> Result<ResolvedPipeline, Report> {
  let text = read_file_to_string(config)?;
  let source = ConfigSource::new(config.display().to_string(), text.clone());

  let value = parse_config_document(&source, &text)?;

  check_pipeline(&source, &value)?;
  let doc = PipelineDoc::from_value(value)?;
  resolve_pipeline(&doc, &process_env())
}

fn process_env() -> Value {
  let env = std::env::vars()
    .map(|(key, value)| (key, Value::String(value)))
    .collect::<Map<_, _>>();
  Value::Object(env)
}

pub(crate) fn run_pipeline(
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

pub(crate) fn select_steps<'a>(
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
        suggestion_suffix(name, &candidates)
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

fn run_step(step: &ResolvedStep, progress: &dyn ProgressSink) -> Result<(), Report> {
  match &step.command {
    PipelineStepCommand::Timetree(args) => {
      let args = TreetimeTimetreeArgs::try_from((**args).clone())?;
      run_timetree_estimation(&args, &NoopCancel, progress).map(|_| ())
    },
    PipelineStepCommand::Optimize(args) => {
      let args = TreetimeOptimizeArgs::try_from(args.clone())?;
      run_optimize(&args, &NoopCancel, progress).map(|_| ())
    },
    PipelineStepCommand::Prune(args) => {
      let args = TreetimePruneArgs::try_from(args.clone())?;
      run_prune(&args, &NoopCancel, progress).map(|_| ())
    },
    PipelineStepCommand::Ancestral(args) => {
      let args = TreetimeAncestralArgs::try_from(args.clone())?;
      run_ancestral_reconstruction(&args, &NoopCancel, progress).map(|_| ())
    },
    PipelineStepCommand::Mugration(args) => {
      let args = TreetimeMugrationArgs::try_from(args.clone())?;
      run_mugration(&args, &NoopCancel, progress).map(|_| ())
    },
    PipelineStepCommand::Clock(args) => {
      let args = TreetimeClockArgs::try_from(args.clone())?;
      run_clock(&args, &NoopCancel, progress).map(|_| ())
    },
  }
}

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

  #[test]
  fn test_runner_select_all_steps_in_order() {
    let pipeline = resolve(config());
    let selected = select_steps(&pipeline, None).unwrap();
    assert_eq!(
      vec!["tt", "anc"],
      selected.iter().map(|step| step.name.as_str()).collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_runner_select_subset_keeps_list_order() {
    let pipeline = resolve(config());
    let selected = select_steps(&pipeline, Some(&btreeset! { "anc".to_owned() })).unwrap();
    assert_eq!(
      vec!["anc"],
      selected.iter().map(|step| step.name.as_str()).collect::<Vec<_>>()
    );
  }

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
